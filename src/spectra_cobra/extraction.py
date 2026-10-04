"""Model extraction.

Both routines work in two phases. First an iterative phase alternates the
forward and reverse LPs over the core reactions, accumulating a flux vector
whose signs say which direction each core reaction should run in. Then a
single network inference problem picks the rest of the network, subject to
those directions.

``spectra_me`` requires a flux consistent model and fails if the core set
cannot all carry flux. ``spectra_ccme`` folds the consistency check into the
same loop, so it accepts an inconsistent model and reports the core reactions
it had to drop.
"""

from logging import getLogger
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Set, Tuple

import numpy as np

from ._lp import carrying_flux, forward, forward_cc, reverse
from ._orientation import (
    STOICHIOMETRY,
    reaction_signs,
    relaxed_mass_balance,
    validate_consistency_type,
)
from .exceptions import SpectraError, SpectraSolverError
from .formulations import (
    MilpSolution,
    growth_optim,
    min_net_lp,
    min_net_milp,
    trade_off,
)

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)

MIN_NET_LP = "minNetLP"
MIN_NET_MILP = "minNetMILP"
GROWTH_OPTIM = "growthOptim"
TRADE_OFF = "tradeOff"
PROBLEM_TYPES = (MIN_NET_LP, MIN_NET_MILP, GROWTH_OPTIM, TRADE_OFF)

#: Problem types that are mixed-integer, and so support pathway exclusion.
MILP_PROBLEM_TYPES = (MIN_NET_MILP, TRADE_OFF)

CORE_DIRECTION = "coreDirection"
PATHWAY_EXCLUSION = "pathwayExclusion"
ALT_SOLUTION_METHODS = (CORE_DIRECTION, PATHWAY_EXCLUSION)

#: The fraction of `tol` at which a core reaction counts as having had its
#: direction settled.
DIRECTION_CUTOFF_FACTOR = 1e-1

#: The range the convex combination weight is drawn from when blending the
#: flux of successive iterations.
BLEND_RANGE = (0.45, 0.55)

#: How far below `tol` a model's solver tolerance should sit. The extraction
#: reads its answer off an LP solution whose mass balance holds only to the
#: solver's tolerance, and a residual of that size breaks the chains carrying
#: flux of order `tol`, which can leave core reactions unable to carry flux.
#: Calibrated on Recon3D: a ratio of 1e3 left dead core reactions, 1e5 left
#: none.
MIN_TOLERANCE_RATIO = 1e4


def _warn_on_loose_tolerance(model: "Model", tol: float) -> None:
    """Warn when the solver's tolerance is too loose for the flux threshold.

    Parameters
    ----------
    model : cobra.Model
        The model about to be extracted from.
    tol : float
        The flux threshold the extraction will use.

    Notes
    -----
    This is the single most effective knob on extraction quality, and it is
    not one most callers would think to touch, so it is worth saying out
    loud rather than leaving in the documentation.

    """
    ratio = tol / model.tolerance if model.tolerance else float("inf")
    if ratio >= MIN_TOLERANCE_RATIO:
        return
    logger.warning(
        "model.tolerance is %.3g against tol=%.3g, a ratio of only %.0f. The "
        "extraction reads its answer off an LP whose mass balance holds to "
        "the solver's tolerance, so a loose one can leave core reactions "
        "present in the result but unable to carry flux. Set "
        "model.tolerance=1e-9, the lowest most solvers accept, before "
        "extracting.",
        model.tolerance,
        tol,
        ratio,
    )


def _normalise_weights(
    model: "Model", weights: Optional[Dict[str, float]]
) -> Dict[str, float]:
    """Return a weight for every reaction, defaulting to one.

    Parameters
    ----------
    model : cobra.Model
        The model whose reactions to weight.
    weights : dict of {str: float}, optional
        The weights to use, keyed by reaction identifier. Reactions left out
        get a weight of 1.0 (default None, i.e. all ones).

    Returns
    -------
    dict of {str: float}
        A weight for each reaction in the model.

    Raises
    ------
    SpectraError
        If `weights` mentions a reaction the model does not have.

    """
    if weights is None:
        return {rxn.id: 1.0 for rxn in model.reactions}

    known = {rxn.id for rxn in model.reactions}
    unknown = set(weights) - known
    if unknown:
        raise SpectraError(
            f"weights refers to {len(unknown)} reaction(s) that are not in the "
            f"model, for example {sorted(unknown)[:5]}."
        )
    return {rxn.id: float(weights.get(rxn.id, 1.0)) for rxn in model.reactions}


def _normalise_core(model: "Model", core_reactions: Iterable) -> Set[str]:
    """Return the core reactions as a set of identifiers.

    Parameters
    ----------
    model : cobra.Model
        The model the core reactions belong to.
    core_reactions : iterable
        The core reactions, given as identifiers, :class:`cobra.Reaction`
        objects, or integer indices into ``model.reactions``.

    Returns
    -------
    set of str
        The identifiers of the core reactions.

    Raises
    ------
    SpectraError
        If any of the given reactions is not in the model.

    """
    core_ids = set()
    for item in core_reactions:
        if isinstance(item, str):
            core_ids.add(item)
        elif isinstance(item, (int, np.integer)):
            core_ids.add(model.reactions[int(item)].id)
        else:
            core_ids.add(item.id)

    known = {rxn.id for rxn in model.reactions}
    unknown = core_ids - known
    if unknown:
        raise SpectraError(
            f"core_reactions refers to {len(unknown)} reaction(s) that are not "
            f"in the model, for example {sorted(unknown)[:5]}."
        )
    return core_ids


def _blend(
    accumulated: Optional[Dict[str, float]],
    new: Dict[str, float],
    rng: np.random.Generator,
) -> Dict[str, float]:
    """Blend a new flux vector into the accumulated one.

    Parameters
    ----------
    accumulated : dict of {str: float}, optional
        The flux accumulated so far, or None on the first iteration.
    new : dict of {str: float}
        The flux just found.
    rng : numpy.random.Generator
        The source of the blending weight.

    Returns
    -------
    dict of {str: float}
        The convex combination ``c * accumulated + (1 - c) * new`` with ``c``
        drawn from [0.45, 0.55], or `new` itself on the first iteration.
        Blending rather than overwriting is what keeps a direction preference
        across iterations.

    """
    if accumulated is None or not any(accumulated.values()):
        return dict(new)
    weight = rng.uniform(*BLEND_RANGE)
    return {
        rxn_id: weight * value + (1.0 - weight) * new[rxn_id]
        for rxn_id, value in accumulated.items()
    }


def _settle_core_directions(
    model: "Model",
    core_ids: Set[str],
    signs: Dict[str, float],
    tol: float,
    rng: np.random.Generator,
    swap_order: bool = False,
) -> Dict[str, float]:
    """Find a flux vector giving every core reaction a direction.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on, already inside a ``with model:`` block and
        with its mass balance relaxed if needed.
    core_ids : set of str
        The core reactions that must end up with a direction.
    signs : dict of {str: float}
        The orientation of each reaction.
    tol : float
        The flux threshold.
    rng : numpy.random.Generator
        The source of random coefficients.
    swap_order : bool, optional
        Whether to try the reverse LP before the forward one, which is how
        the ``coreDirection`` alternative-solution method explores different
        directions (default False).

    Returns
    -------
    dict of {str: float}
        The accumulated flux vector.

    Raises
    ------
    SpectraInfeasibleCoreError
        If the core reactions cannot all carry flux at once.

    """
    unsettled = set(core_ids)
    accumulated = None
    cutoff = tol * DIRECTION_CUTOFF_FACTOR
    first, second = (reverse, forward) if swap_order else (forward, reverse)

    while unsettled:
        with model:
            fluxes = first(model, unsettled, signs, tol, rng)
        if fluxes is None:
            # Only reverse() can decline to raise; treat it as no progress.
            fluxes = {rxn.id: 0.0 for rxn in model.reactions}
        accumulated = _blend(accumulated, fluxes, rng)
        unsettled -= carrying_flux({r: accumulated[r] for r in unsettled}, cutoff)
        if not unsettled:
            break

        with model:
            fluxes = second(model, unsettled, signs, tol, rng)
        if fluxes is None:
            fluxes = {rxn.id: 0.0 for rxn in model.reactions}
        accumulated = _blend(accumulated, fluxes, rng)
        settled = carrying_flux({r: accumulated[r] for r in unsettled}, cutoff)
        if not settled:
            raise SpectraError(
                f"{len(unsettled)} core reaction(s) could not be given a flux "
                f"direction, for example {sorted(unsettled)[:5]}. Their "
                f"accumulated flux stays below {cutoff:.3g}, which can happen "
                f"when the convex combination of two iterations cancels out; "
                f"retry with a different seed."
            )
        unsettled -= settled

    return accumulated


def _directions_from_flux(
    model: "Model",
    fluxes: Dict[str, float],
    core_ids: Set[str],
    signs: Dict[str, float],
) -> Dict[str, int]:
    """Read the oriented direction of each core reaction off a flux vector.

    Parameters
    ----------
    model : cobra.Model
        The model the fluxes belong to.
    fluxes : dict of {str: float}
        The accumulated flux vector.
    core_ids : set of str
        The reactions whose direction matters.
    signs : dict of {str: float}
        The orientation of each reaction.

    Returns
    -------
    dict of {str: int}
        1 where the oriented flux is positive, -1 where it is negative, and 0
        elsewhere.

    """
    directions = {}
    for rxn in model.reactions:
        if rxn.id not in core_ids:
            directions[rxn.id] = 0
            continue
        oriented = signs[rxn.id] * fluxes.get(rxn.id, 0.0)
        directions[rxn.id] = 1 if oriented > 0 else (-1 if oriented < 0 else 0)
    return directions


def _solve_formulation(
    model: "Model",
    directions: Dict[str, int],
    weights: Dict[str, float],
    tol: float,
    steady_state: bool,
    problem_type: str,
    time_limit: Optional[float],
    previous_solutions: Optional[List[Set[str]]] = None,
    inclusion_cutoff: Optional[float] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
) -> MilpSolution:
    """Dispatch to the requested network inference formulation.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from.
    directions : dict of {str: int}
        The oriented direction each reaction must carry flux in.
    weights : dict of {str: float}
        The weight of each reaction.
    tol : float
        The flux threshold.
    steady_state : bool
        Whether to enforce ``S v = 0`` rather than ``S v >= 0``.
    problem_type : str
        One of :data:`PROBLEM_TYPES`.
    time_limit : float, optional
        The time limit for the mixed-integer formulations.
    previous_solutions : list of set of str, optional
        Reaction sets to exclude, for the mixed-integer formulations
        (default None).

    Returns
    -------
    MilpSolution
        The reactions to keep, and the ones a mixed-integer formulation's
        indicators selected. The linear formulations have no indicators, so
        the two are the same for them.

    Raises
    ------
    SpectraError
        If `problem_type` is unknown, or if `previous_solutions` is given for
        a formulation that cannot honour it.

    """
    if problem_type not in PROBLEM_TYPES:
        raise SpectraError(
            f"problem_type must be one of {PROBLEM_TYPES}, not {problem_type!r}."
        )
    if previous_solutions and problem_type not in MILP_PROBLEM_TYPES:
        raise SpectraError(
            f"Excluding previous solutions needs a mixed-integer formulation, "
            f"one of {MILP_PROBLEM_TYPES}, not {problem_type!r}."
        )

    if indicator_reactions is not None and problem_type not in MILP_PROBLEM_TYPES:
        raise SpectraError(
            f"indicator_reactions needs a mixed-integer formulation, one of "
            f"{MILP_PROBLEM_TYPES}, not {problem_type!r}."
        )

    if problem_type == MIN_NET_LP:
        kept = min_net_lp(
            model, directions, weights, tol, steady_state, inclusion_cutoff
        )
        return MilpSolution(included=kept, selected=kept)
    if problem_type == GROWTH_OPTIM:
        kept = growth_optim(
            model, directions, weights, tol, steady_state, inclusion_cutoff
        )
        return MilpSolution(included=kept, selected=kept)
    if problem_type == MIN_NET_MILP:
        return min_net_milp(
            model,
            directions,
            weights,
            tol,
            steady_state,
            time_limit,
            previous_solutions,
            indicator_reactions,
            inclusion_cutoff,
        )
    return trade_off(
        model,
        directions,
        weights,
        tol,
        steady_state,
        time_limit,
        previous_solutions,
        indicator_reactions,
        inclusion_cutoff,
    )


def _extract(model: "Model", keep_ids: Set[str], remove_genes: bool) -> "Model":
    """Return a copy of the model holding only the given reactions.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from.
    keep_ids : set of str
        The identifiers of the reactions to keep.
    remove_genes : bool
        Whether to drop the genes left without a reaction.

    Returns
    -------
    cobra.Model
        The extracted model.

    Notes
    -----
    Orphaned metabolites are always dropped, orphaned genes only when asked.
    cobrapy's ``remove_orphans`` covers both at once, so it is left off here
    and the two are handled separately; otherwise ``remove_genes=False``
    would remove the genes regardless.

    """
    extracted = model.copy()
    extracted.remove_reactions(
        [rxn.id for rxn in model.reactions if rxn.id not in keep_ids],
        remove_orphans=False,
    )
    orphaned_metabolites = [met for met in extracted.metabolites if not met.reactions]
    if orphaned_metabolites:
        extracted.remove_metabolites(orphaned_metabolites)
    if remove_genes:
        from cobra.manipulation import remove_genes as _remove_genes

        unused = [gene.id for gene in extracted.genes if not gene.reactions]
        if unused:
            _remove_genes(extracted, unused, remove_reactions=False)
    return extracted


def spectra_me(
    model: "Model",
    core_reactions: Iterable,
    tol: float = 1e-4,
    consistency_type: str = STOICHIOMETRY,
    weights: Optional[Dict[str, float]] = None,
    n_solutions: int = 1,
    alt_solution_method: str = CORE_DIRECTION,
    problem_type: str = MIN_NET_LP,
    time_limit: Optional[float] = 7200.0,
    remove_genes: bool = True,
    previous_solutions: Optional[List[Set[str]]] = None,
    seed: Optional[int] = None,
    inclusion_cutoff: Optional[float] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
) -> "Model":
    """Extract a context-specific model around a set of core reactions.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from. It must be flux consistent; run
        :func:`spectra_cc` first, or use :func:`spectra_ccme` instead.
    core_reactions : iterable
        The reactions that must be in the extracted model, as identifiers,
        :class:`cobra.Reaction` objects, or indices.
    tol : float, optional
        The minimum absolute flux every reaction in the extracted model has to
        carry (default 1e-4).
    consistency_type : {"stoichiometry", "topology"}, optional
        Whether to assume a steady state or an accumulation condition
        (default "stoichiometry").
    weights : dict of {str: float}, optional
        The weight of each reaction, keyed by identifier; reactions left out
        get 1.0. What a weight means depends on `problem_type`: the
        ``minNet`` formulations and ``growthOptim`` treat it as a cost, so it
        should be non-negative, while ``tradeOff`` treats it as evidence and
        accepts any real number (default None, i.e. all ones).
    n_solutions : int, optional
        How many alternative models to return (default 1).
    alt_solution_method : {"coreDirection", "pathwayExclusion"}, optional
        How to find alternatives. ``"coreDirection"`` re-runs the iterative
        phase with the two LPs in a random order, giving the core reactions
        different directions; ``"pathwayExclusion"`` forbids each solution
        already found and so needs a mixed-integer `problem_type`
        (default "coreDirection").
    problem_type : str, optional
        The network inference formulation, one of :data:`PROBLEM_TYPES`
        (default "minNetLP").
    time_limit : float, optional
        The maximum time to spend on each mixed-integer solve, in seconds
        (default 7200).
    inclusion_cutoff : float, optional
        The absolute flux at which a reaction counts as part of the extracted
        model, for the LP formulations (default ``tol * 1e-7``). The default
        is far below the solver's own tolerance on purpose: raising it drops
        reactions that carry almost no flux but are load-bearing for a mass
        balance, which leaves core reactions present in the result yet unable
        to carry flux. The mixed-integer formulations ignore it, since they
        read their answer off their binaries.
    remove_genes : bool, optional
        Whether to drop the genes left without a reaction by the extraction
        (default True). Orphaned metabolites are always dropped.
    previous_solutions : list of set of str, optional
        Reaction sets to exclude from the first solution too, for a
        mixed-integer `problem_type` (default None).
    seed : int, optional
        A seed for the random coefficients, making the result reproducible
        (default None).

    Returns
    -------
    cobra.Model or list of cobra.Model
        The extracted model, or a list of `n_solutions` models if more than
        one was asked for.

    Raises
    ------
    SpectraInfeasibleCoreError
        If the core reactions cannot all carry at least `tol` of flux at once,
        which is what happens when any of them is blocked in `model`.
    SpectraError
        If an argument is invalid, or if a core reaction cannot be given a
        direction.

    See Also
    --------
    spectra_ccme : The same, for a model that is not known to be consistent.

    """
    _warn_on_loose_tolerance(model, tol)
    models, _ = _spectra_me(
        model,
        core_reactions,
        tol=tol,
        consistency_type=consistency_type,
        weights=weights,
        n_solutions=n_solutions,
        alt_solution_method=alt_solution_method,
        problem_type=problem_type,
        time_limit=time_limit,
        remove_genes=remove_genes,
        previous_solutions=previous_solutions,
        seed=seed,
        inclusion_cutoff=inclusion_cutoff,
        indicator_reactions=indicator_reactions,
    )
    return models[0] if n_solutions == 1 else models


def _spectra_me(
    model: "Model",
    core_reactions: Iterable,
    tol: float,
    consistency_type: str,
    weights: Optional[Dict[str, float]],
    n_solutions: int,
    alt_solution_method: str,
    problem_type: str,
    time_limit: Optional[float],
    remove_genes: bool,
    previous_solutions: Optional[List[Set[str]]],
    seed: Optional[int],
    blocked_ids: Optional[Set[str]] = None,
    inclusion_cutoff: Optional[float] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
) -> Tuple[List["Model"], List[Set[str]]]:
    """Run the extraction, returning the models and the reaction sets found.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from.
    core_reactions : iterable
        The core reactions.
    tol : float
        The flux threshold.
    consistency_type : str
        Whether to assume a steady state or an accumulation condition.
    weights : dict of {str: float}, optional
        The weight of each reaction.
    n_solutions : int
        How many models to return.
    alt_solution_method : str
        How to find alternatives.
    problem_type : str
        The network inference formulation.
    time_limit : float, optional
        The time limit for mixed-integer solves.
    remove_genes : bool
        Whether to drop orphaned genes.
    previous_solutions : list of set of str, optional
        Reaction sets to exclude.
    seed : int, optional
        A seed for the random coefficients.
    blocked_ids : set of str, optional
        Core reactions already known to be blocked, which are dropped from the
        core set rather than being required to carry flux (default None).

    Returns
    -------
    tuple of (list of cobra.Model, list of set of str)
        The extracted models and the reaction sets they were built from.

    """
    if n_solutions < 1:
        raise SpectraError(f"n_solutions must be at least 1, not {n_solutions}.")
    if alt_solution_method not in ALT_SOLUTION_METHODS:
        raise SpectraError(
            f"alt_solution_method must be one of {ALT_SOLUTION_METHODS}, "
            f"not {alt_solution_method!r}."
        )

    steady_state = validate_consistency_type(consistency_type)
    all_weights = _normalise_weights(model, weights)
    core_ids = _normalise_core(model, core_reactions)
    if blocked_ids:
        core_ids -= blocked_ids
    signs = reaction_signs(model)
    rng = np.random.default_rng(seed)

    exclude = [set(s) for s in previous_solutions] if previous_solutions else []
    models: List["Model"] = []
    found: List[Set[str]] = []

    for index in range(n_solutions):
        # The first solution always runs the LPs in their canonical order; the
        # alternatives shuffle it, which is what makes them differ.
        swap_order = (
            index > 0
            and alt_solution_method == CORE_DIRECTION
            and bool(rng.integers(2))
        )

        if index == 0 or alt_solution_method == CORE_DIRECTION:
            with model, relaxed_mass_balance(model, steady_state):
                fluxes = _settle_core_directions(
                    model, core_ids, signs, tol, rng, swap_order
                )
            directions = _directions_from_flux(model, fluxes or {}, core_ids, signs)

        try:
            solution = _solve_formulation(
                model,
                directions,
                all_weights,
                tol,
                steady_state,
                problem_type,
                time_limit,
                exclude or None,
                inclusion_cutoff,
                indicator_reactions,
            )
        except SpectraSolverError:
            # Excluding every solution found so far can leave the problem with
            # no feasible answer, which simply means the network has no further
            # alternative to offer, so the solutions already found still
            # stand. A failure on the very first solve is a real error, so it
            # is left to propagate.
            if index > 0 and alt_solution_method == PATHWAY_EXCLUSION:
                logger.info(
                    "No further alternative solution exists; returning the %d "
                    "found rather than the %d requested.",
                    index,
                    n_solutions,
                )
                break
            raise

        models.append(_extract(model, solution.included, remove_genes))
        found.append(solution.included)

        if alt_solution_method == PATHWAY_EXCLUSION:
            # The exclusion constraint is written over the indicators, so it
            # has to be given the reactions they actually selected. Handing
            # it the wider kept set would let the same solution recur: the
            # extra members sit at zero, so the sum stays under the bound.
            exclude.append(solution.selected)

    return models, found


def spectra_ccme(
    model: "Model",
    core_reactions: Iterable,
    tol: float = 1e-4,
    consistency_type: str = STOICHIOMETRY,
    weights: Optional[Dict[str, float]] = None,
    n_solutions: int = 1,
    alt_solution_method: str = CORE_DIRECTION,
    problem_type: str = MIN_NET_LP,
    time_limit: Optional[float] = 7200.0,
    remove_genes: bool = True,
    seed: Optional[int] = None,
    inclusion_cutoff: Optional[float] = None,
) -> Tuple["Model", List[str]]:
    """Check consistency and extract a model in a single pass.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from. It need not be flux consistent.
    core_reactions : iterable
        The reactions that must be in the extracted model, as identifiers,
        :class:`cobra.Reaction` objects, or indices. Any of them that is
        blocked is dropped and reported rather than causing a failure.
    tol : float, optional
        The minimum absolute flux required for a reaction to count as
        unblocked, and to be carried by the extracted model (default 1e-4).
    consistency_type : {"stoichiometry", "topology"}, optional
        Whether to assume a steady state or an accumulation condition
        (default "stoichiometry").
    weights : dict of {str: float}, optional
        The weight of each reaction (default None, i.e. all ones).
    n_solutions : int, optional
        How many alternative models to return (default 1).
    alt_solution_method : {"coreDirection", "pathwayExclusion"}, optional
        How to find alternatives (default "coreDirection").
    problem_type : str, optional
        The network inference formulation, one of :data:`PROBLEM_TYPES`
        (default "minNetLP").
    time_limit : float, optional
        The maximum time to spend on each mixed-integer solve, in seconds
        (default 7200).
    inclusion_cutoff : float, optional
        The absolute flux at which a reaction counts as part of the extracted
        model, for the LP formulations (default ``tol * 1e-7``). The default
        is far below the solver's own tolerance on purpose: raising it drops
        reactions that carry almost no flux but are load-bearing for a mass
        balance, which leaves core reactions present in the result yet unable
        to carry flux. The mixed-integer formulations ignore it, since they
        read their answer off their binaries.
    remove_genes : bool, optional
        Whether to drop the genes left without a reaction by the extraction
        (default True). Orphaned metabolites are always dropped.
    seed : int, optional
        A seed for the random coefficients (default None).

    Returns
    -------
    tuple of (cobra.Model or list of cobra.Model, list of str)
        The extracted model (or a list of them if `n_solutions` is above one)
        and the identifiers of the core reactions that turned out to be
        blocked and so are absent from it.

    See Also
    --------
    spectra_cc : The consistency check on its own.
    spectra_me : Extraction from a model already known to be consistent.

    """
    steady_state = validate_consistency_type(consistency_type)
    _warn_on_loose_tolerance(model, tol)
    core_ids = _normalise_core(model, core_reactions)
    signs = reaction_signs(model)
    rng = np.random.default_rng(seed)

    # The consistency loop of spectra_cc, but keeping the flux it finds so the
    # core directions can be read off it rather than recomputed.
    rxns_to_check = {rxn.id for rxn in model.reactions}
    accumulated = None
    cutoff = tol * 0.99

    with model, relaxed_mass_balance(model, steady_state):
        previous_count = None
        while len(rxns_to_check) != previous_count:
            previous_count = len(rxns_to_check)

            with model:
                fluxes = forward_cc(model, rxns_to_check, signs, tol, rng)
            if fluxes is not None:
                rxns_to_check -= carrying_flux(fluxes, cutoff)
                accumulated = _blend(accumulated, fluxes, rng)

            with model:
                fluxes = reverse(model, rxns_to_check, signs, tol, rng)
            if fluxes is not None:
                rxns_to_check -= carrying_flux(fluxes, cutoff)
                accumulated = _blend(accumulated, fluxes, rng)

    blocked_ids = set(rxns_to_check)
    blocked_core = [
        rxn.id for rxn in model.reactions if rxn.id in core_ids & blocked_ids
    ]
    if blocked_core:
        logger.warning(
            "%d of the %d core reactions are blocked and were dropped: %s",
            len(blocked_core),
            len(core_ids),
            blocked_core[:5],
        )

    live_core = core_ids - blocked_ids
    directions = _directions_from_flux(model, accumulated or {}, live_core, signs)
    undirected = [r for r in live_core if directions[r] == 0]
    if undirected:
        raise SpectraError(
            f"{len(undirected)} unblocked core reaction(s) ended up with zero "
            f"accumulated flux, for example {sorted(undirected)[:5]}, so their "
            f"direction is unknown. This happens when the convex combination "
            f"of two iterations cancels out; retry with a different seed."
        )

    all_weights = _normalise_weights(model, weights)
    solution = _solve_formulation(
        model,
        directions,
        all_weights,
        tol,
        steady_state,
        problem_type,
        time_limit,
        None,
        inclusion_cutoff,
    )
    first = _extract(model, solution.included, remove_genes)

    if n_solutions == 1:
        return first, blocked_core

    # Further solutions come from spectra_me on the consistent sub-model,
    # excluding the solution just found.
    consistent = model.copy()
    consistent.remove_reactions(sorted(blocked_ids), remove_orphans=True)
    consistent_ids = {rxn.id for rxn in consistent.reactions}
    rest, _ = _spectra_me(
        consistent,
        sorted(live_core),
        tol=tol,
        consistency_type=consistency_type,
        weights={k: v for k, v in all_weights.items() if k in consistent_ids},
        n_solutions=n_solutions - 1,
        alt_solution_method=alt_solution_method,
        problem_type=problem_type,
        time_limit=time_limit,
        remove_genes=remove_genes,
        previous_solutions=(
            # The indicator-selected set, not the wider kept set: see the
            # note where _spectra_me builds its own exclusion list.
            [solution.selected & consistent_ids]
            if alt_solution_method == PATHWAY_EXCLUSION
            else None
        ),
        seed=None if seed is None else seed + 1,
    )
    return [first] + rest, blocked_core
