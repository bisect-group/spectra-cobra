"""Model extraction (`spectra_me`, `spectra_ccme`, and `spectra_me2`).

All extraction routines operate in two stages without any convex combinations
of LP flux vectors:

1. **Stage 1 (Core Direction Discovery):** Classifies core reactions in 0 LPs
   by bound reachability (`can_fwd`, `can_rev`), pins one-way core reactions
   immediately, runs a 1-LP Split-Variable `minNetLP` Anchor when fixed flux
   exists (including `-weights` penalty steering for negative-weight `tradeOff`),
   and settles any remaining two-way core reactions via iterative `mPinned`
   `forward`/`reverse` LPs sharing per-round random weights.
2. **Stage 2 (Network Inference & Certification):** Dispatches to the selected
   formulation (`minNetLP`, `minNetMILP`, `tradeOff`, or `optimBiomass`) with
   0-solve `zero_feas` skips and support consistency certification.
"""

import math
import warnings
from logging import getLogger
from typing import (
    TYPE_CHECKING,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Set,
    Tuple,
    Union,
)

import numpy as np
from cobra.util.solver import linear_reaction_coefficients
from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

from ._copy import subset_model
from ._lp import WEIGHT_RANGE, forward, forward_cc, reverse
from ._orientation import (
    STOICHIOMETRY,
    apply_direction_bounds,
    oriented_bounds,
    reaction_signs,
    relaxed_mass_balance,
    validate_consistency_type,
)
from .exceptions import (
    SpectraError,
    SpectraInfeasibleCoreError,
    SpectraSolverError,
    SpectraValueError,
)
from .formulations import (
    DEFAULT_SOLVER_TOLERANCE,
    MilpSolution,
    _LpToleranceConfig,
    _apply_tolerances,
    _build_network_arrays,
    _certify_milp_support,
    _min_net_lp_core,
    _min_net_milp_core,
    _optim_biomass_core,
    _resolve_lp_tolerances,
    _strip_dead_ends,
    _trade_off_core,
    _validate_positive_finite,
    _zero_flux_is_feasible,
)

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)

MIN_NET_LP = "minNetLP"
MIN_NET_MILP = "minNetMILP"
OPTIM_BIOMASS = "optimBiomass"
TRADE_OFF = "tradeOff"
PROBLEM_TYPES = (MIN_NET_LP, MIN_NET_MILP, OPTIM_BIOMASS, TRADE_OFF)

#: Problem types that are mixed-integer, and so support pathway exclusion.
MILP_PROBLEM_TYPES = (MIN_NET_MILP, TRADE_OFF)

CORE_DIRECTION = "coreDirection"
PATHWAY_EXCLUSION = "pathwayExclusion"
ALT_SOLUTION_METHODS = (CORE_DIRECTION, PATHWAY_EXCLUSION)

#: Fraction of `tol` at which the 1-LP Stage-1 Anchor tentatively settles a
#: two-way core reaction's sign.
DIRECTION_CUTOFF_FACTOR = 1e-1

#: Maximum consecutive stalled rounds in the iterative direction-settling loop.
DIRECTION_ATTEMPTS = 5

#: Minimum recommended ratio `tol / feasibility_tolerance`.
MIN_TOLERANCE_RATIO = 1e4


def _warn_on_loose_tolerance(
    model: "Model", tol: float, feasibility: Optional[float] = None
) -> None:
    """Warn when the solver's tolerance is too loose for the flux threshold."""
    effective = model.tolerance if feasibility is None else feasibility
    ratio = tol / effective if effective else float("inf")
    if ratio >= MIN_TOLERANCE_RATIO:
        return
    logger.warning(
        "The solver's feasibility tolerance is %.3g against tol=%.3g, a ratio "
        "of only %.0f. The extraction reads its answer off an LP whose mass "
        "balance holds to that tolerance, so a loose one can leave core "
        "reactions present in the result but unable to carry flux. Pass "
        "feasibility_tolerance=1e-9, the lowest most solvers accept.",
        effective,
        tol,
        ratio,
    )


def _validate_problem_type(problem_type: str) -> str:
    """Validate `problem_type` against the four supported v2 formulations."""
    if problem_type not in PROBLEM_TYPES:
        raise SpectraValueError(
            f"problem_type must be one of {PROBLEM_TYPES}, not {problem_type!r}."
        )
    return problem_type


def _resolve_alt_solution_method(
    alt_solution_method: Optional[str],
    problem_type: str,
    n_solutions: int,
) -> str:
    """Resolve default `alt_solution_method` by formulation class and validate if `n_solutions > 1`."""
    is_milp = problem_type in MILP_PROBLEM_TYPES
    if alt_solution_method is None:
        resolved = PATHWAY_EXCLUSION if is_milp else CORE_DIRECTION
    else:
        resolved = alt_solution_method

    if n_solutions > 1:
        if resolved not in ALT_SOLUTION_METHODS:
            raise SpectraValueError(
                f"alt_solution_method must be one of {ALT_SOLUTION_METHODS}, "
                f"not {resolved!r}."
            )
        if resolved == PATHWAY_EXCLUSION and not is_milp:
            raise SpectraValueError(
                f"alt_solution_method {PATHWAY_EXCLUSION!r} needs a mixed-integer "
                f"formulation, one of {MILP_PROBLEM_TYPES}, not {problem_type!r}."
            )
    return resolved


def _normalise_weights(
    model: "Model",
    weights: Optional[Mapping[str, float]],
    problem_type: str = MIN_NET_LP,
    indicator_groups: Optional[Mapping[str, Iterable[str]]] = None,
) -> Dict[str, float]:
    """Return validated weights for reactions and optional indicator groups."""
    allow_negative = problem_type == TRADE_OFF
    rxn_ids = [rxn.id for rxn in model.reactions]
    rxn_set = set(rxn_ids)
    extra = set(indicator_groups or ())

    if weights is None:
        all_weights = {rid: 1.0 for rid in rxn_ids}
        all_weights.update({key: 1.0 for key in extra})
        return all_weights

    known = rxn_set | extra
    unknown = set(weights) - known
    if unknown:
        raise SpectraValueError(
            f"weights refers to {len(unknown)} reaction(s) that are not in the "
            f"model, for example {sorted(unknown)[:5]}."
        )

    for key, raw in weights.items():
        if isinstance(raw, bool) or not isinstance(raw, (int, float, np.number)):
            raise SpectraValueError(
                f"Weight for {key!r} must be a finite numeric value, got {raw!r}."
            )
        val = float(raw)
        if not math.isfinite(val):
            raise SpectraValueError(
                f"Weight for {key!r} must be finite, got {raw!r}."
            )
        if not allow_negative and val < 0.0:
            raise SpectraValueError(
                f"weights must be non-negative for problem_type {problem_type!r} "
                f"(got {val:g} on {key!r}); use 'tradeOff' for signed weights."
            )

    all_weights = {rid: float(weights.get(rid, 1.0)) for rid in rxn_ids}
    if indicator_groups is not None:
        for key, members in indicator_groups.items():
            if key in weights:
                all_weights[key] = float(weights[key])
            else:
                member_list = list(members)
                if member_list and any(m in weights for m in member_list):
                    all_weights[key] = sum(
                        float(weights.get(m, 1.0)) for m in member_list
                    )
                else:
                    all_weights[key] = 1.0
    return all_weights


def _normalise_core(model: "Model", core_reactions: Iterable) -> Set[str]:
    """Return the core reactions as a validated set of identifiers."""
    core_ids: Set[str] = set()
    n_rxns = len(model.reactions)
    for item in core_reactions:
        if isinstance(item, str):
            core_ids.add(item)
        elif isinstance(item, (int, np.integer)) and not isinstance(item, bool):
            idx = int(item)
            if idx < 0 or idx >= n_rxns:
                raise SpectraValueError(
                    f"Core reaction index {idx} is out of range for model with {n_rxns} reactions."
                )
            core_ids.add(model.reactions[idx].id)
        else:
            core_ids.add(item.id)

    known = {rxn.id for rxn in model.reactions}
    unknown = core_ids - known
    if unknown:
        raise SpectraValueError(
            f"core_reactions refers to {len(unknown)} reaction(s) that are not "
            f"in the model, for example {sorted(unknown)[:5]}."
        )
    return core_ids


def _raise_infeasible_core(n_core: int, tol: float) -> None:
    """Raise `SpectraInfeasibleCoreError` with actionable guidance."""
    raise SpectraInfeasibleCoreError(
        f"No flux distribution lets all {n_core} core reactions carry at least "
        f"tol={tol:.3g} at once. Either the core set is mutually inconsistent, "
        f"or some of its reactions are blocked in this model; spectra_ccme "
        f"handles the latter by dropping them and returning them."
    )


def _assert_core_progress(
    unsettled: Set[str],
    before: int,
    stalls: int,
    tol: float,
    core_is_known_live: bool = False,
) -> int:
    """Track stalled rounds and raise if `DIRECTION_ATTEMPTS` is exceeded."""
    if len(unsettled) < before:
        return 0
    stalls += 1
    if stalls > DIRECTION_ATTEMPTS:
        detail = (
            ""
            if core_is_known_live
            else (
                " Either they are blocked in the model as given -- "
                "spectra_me requires a flux consistent model, so run "
                "spectra_cc first or use spectra_ccme, which drops a "
                "blocked core reaction instead of failing -- or the "
                "weights were unlucky, in which case another seed settles "
                "them."
            )
        )
        raise SpectraError(
            f"{len(unsettled)} core reaction(s) could not be given a flux "
            f"direction in {DIRECTION_ATTEMPTS} attempts, for example "
            f"{sorted(unsettled)[:5]} at tol={tol:.3g}. Reactions that can "
            f"only carry flux together, in opposite directions, are the "
            f"usual cause.{detail}"
        )
    return stalls


def _discover_core_directions(
    model: "Model",
    core_ids: Set[str],
    signs: Dict[str, float],
    all_weights: Mapping[str, float],
    problem_type: str,
    tol: float,
    steady_state: bool,
    inclusion_cutoff: Optional[float],
    feasibility_tolerance: Optional[float],
    optimality_tolerance: Optional[float],
    lock_tol: float,
    rng: np.random.Generator,
    is_alt_sol: bool = False,
    has_per_rxn_weights: bool = True,
    core_is_known_live: bool = False,
) -> Tuple[Dict[str, int], Optional[Dict[str, float]], Optional[Dict[str, float]], int]:
    """Discover simultaneously feasible core reaction directions with 0 convex combinations."""
    directions = {rxn.id: 0 for rxn in model.reactions}
    last_flux: Optional[Dict[str, float]] = None
    pass0_flux_ready: Optional[Dict[str, float]] = None
    n_lps = 0

    # Skip 1: Empty core set -> 0 LPs
    if not core_ids:
        return directions, last_flux, pass0_flux_ready, n_lps

    # Classify core reactions by bound reachability at lock_tol in 0 LPs
    only_fwd_core: Set[str] = set()
    only_rev_core: Set[str] = set()
    two_way_core: Set[str] = set()
    for rxn_id in core_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        lo, hi = oriented_bounds(rxn, signs[rxn_id])
        can_fwd = hi >= lock_tol
        can_rev = lo <= -lock_tol
        if not can_fwd and not can_rev:
            _raise_infeasible_core(len(core_ids), tol)
        if can_fwd and not can_rev:
            only_fwd_core.add(rxn_id)
            directions[rxn_id] = 1
        elif can_rev and not can_fwd:
            only_rev_core.add(rxn_id)
            directions[rxn_id] = -1
        else:
            two_way_core.add(rxn_id)

    # Skip 2: All core reactions are one-way -> 0 Stage-1 LPs
    if not two_way_core:
        return directions, last_flux, pass0_flux_ready, n_lps

    # Build non-negative Stage-1 Anchor penalty vector when fixed flux exists
    has_fixed = bool(
        only_fwd_core
        or only_rev_core
        or not _zero_flux_is_feasible(model, steady_state)
    )
    use_anchor = False
    w_anchor: Optional[Dict[str, float]] = None
    if (not is_alt_sol) and has_per_rxn_weights and has_fixed:
        if problem_type in (MIN_NET_LP, MIN_NET_MILP):
            w_anchor = {
                rxn.id: (0.0 if rxn.id in core_ids else float(all_weights[rxn.id]))
                for rxn in model.reactions
            }
            use_anchor = True
        elif problem_type == TRADE_OFF and any(
            float(all_weights[rxn.id]) < 0.0 for rxn in model.reactions
        ):
            w_anchor = {}
            for rxn in model.reactions:
                if rxn.id in core_ids:
                    w_anchor[rxn.id] = 0.0
                else:
                    w_val = float(all_weights[rxn.id])
                    w_anchor[rxn.id] = -w_val if w_val < 0.0 else 0.0
            use_anchor = True

    v_anchor: Dict[str, float] = {rxn.id: 0.0 for rxn in model.reactions}
    if use_anchor and w_anchor is not None:
        try:
            _, v_anchor, anc_lps = _min_net_lp_core(
                model=model,
                directions=directions,
                weights=w_anchor,
                tol=tol,
                steady_state=steady_state,
                inclusion_cutoff=inclusion_cutoff,
                feasibility_tolerance=feasibility_tolerance,
                optimality_tolerance=optimality_tolerance,
                run_irls=False,
            )
        except SpectraSolverError as err:
            raise SpectraInfeasibleCoreError(
                f"No flux distribution lets all {len(core_ids)} core reactions "
                f"carry at least tol={tol:.3g} at once. Use spectra_ccme if "
                f"the input model may contain blocked reactions."
            ) from err

        n_lps += anc_lps
        last_flux = v_anchor

        for rxn_id in two_way_core:
            if abs(v_anchor[rxn_id]) >= DIRECTION_CUTOFF_FACTOR * tol:
                ori = signs[rxn_id] * v_anchor[rxn_id]
                directions[rxn_id] = 1 if ori > 0.0 else -1

        unsettled = {r for r in two_way_core if directions[r] == 0}
        if not unsettled:
            if all(abs(v_anchor[r]) >= tol for r in core_ids):
                if problem_type == MIN_NET_LP:
                    pass0_flux_ready = v_anchor
                return directions, last_flux, pass0_flux_ready, n_lps
            if all(abs(v_anchor[r]) >= lock_tol for r in core_ids):
                return directions, last_flux, pass0_flux_ready, n_lps
            for rxn_id in two_way_core:
                if abs(v_anchor[rxn_id]) < lock_tol:
                    directions[rxn_id] = 0
            unsettled = {r for r in two_way_core if directions[r] == 0}
    else:
        unsettled = set(two_way_core)

    # Iterative mPinned forward/reverse on unsettled two-way core reactions
    stalls = 0
    round_idx = 0
    with model, relaxed_mass_balance(model, steady_state):
        while unsettled:
            round_idx += 1
            before = len(unsettled)
            first_is_fwd = bool(rng.integers(2) == 0) if is_alt_sol else True

            ordered_u = sorted(unsettled)
            iter_weights = dict(
                zip(ordered_u, rng.uniform(*WEIGHT_RANGE, size=len(ordered_u)))
            )

            pass_order = (first_is_fwd, not first_is_fwd)
            for p_idx, do_fwd in enumerate(pass_order):
                if not unsettled:
                    break

                def _solve_pinned() -> Optional[Dict[str, float]]:
                    with model:
                        apply_direction_bounds(model, directions, signs, tol)
                        if do_fwd:
                            return forward(
                                model,
                                unsettled,
                                signs,
                                tol,
                                rng,
                                objective_weights=iter_weights,
                            )
                        return reverse(
                            model,
                            unsettled,
                            signs,
                            tol,
                            rng,
                            raise_on_failure=True,
                            objective_weights=iter_weights,
                        )

                n_lps += 1
                try:
                    f = _solve_pinned()
                except (SpectraInfeasibleCoreError, SpectraSolverError) as err:
                    if use_anchor and round_idx == 1 and p_idx == 0:
                        for rxn_id in two_way_core:
                            if abs(v_anchor[rxn_id]) < lock_tol:
                                directions[rxn_id] = 0
                        unsettled = {r for r in two_way_core if directions[r] == 0}
                        ordered_u = sorted(unsettled)
                        iter_weights = dict(
                            zip(
                                ordered_u,
                                rng.uniform(*WEIGHT_RANGE, size=len(ordered_u)),
                            )
                        )
                        n_lps += 1
                        try:
                            f = _solve_pinned()
                        except (SpectraInfeasibleCoreError, SpectraSolverError) as err2:
                            raise SpectraInfeasibleCoreError(
                                f"No flux distribution lets all {len(core_ids)} core "
                                f"reactions carry at least tol={tol:.3g} at once. "
                                f"Use spectra_ccme if the model may contain blocked reactions."
                            ) from err2
                    else:
                        raise SpectraInfeasibleCoreError(
                            f"No flux distribution lets all {len(core_ids)} core "
                            f"reactions carry at least tol={tol:.3g} at once. "
                            f"Use spectra_ccme if the model may contain blocked reactions."
                        ) from err

                if f is None:
                    _raise_infeasible_core(len(core_ids), tol)

                last_flux = f
                for rxn_id in list(unsettled):
                    if abs(f[rxn_id]) >= lock_tol:
                        ori = signs[rxn_id] * f[rxn_id]
                        directions[rxn_id] = 1 if ori > 0.0 else -1
                        unsettled.remove(rxn_id)

            stalls = _assert_core_progress(
                unsettled, before, stalls, tol, core_is_known_live=core_is_known_live
            )

    return directions, last_flux, pass0_flux_ready, n_lps


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
    indicator_groups: Optional[Mapping[str, Iterable[str]]] = None,
    feasibility_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    optimality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    integrality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    warm_start: Optional[Mapping[str, float]] = None,
    pass0_flux_ready: Optional[Mapping[str, float]] = None,
) -> Tuple[MilpSolution, Optional[np.ndarray], int]:
    """Dispatch to the requested network inference formulation."""
    _validate_problem_type(problem_type)

    if previous_solutions and problem_type not in MILP_PROBLEM_TYPES:
        raise SpectraValueError(
            f"Excluding previous solutions needs a mixed-integer formulation, "
            f"one of {MILP_PROBLEM_TYPES}, not {problem_type!r}."
        )
    if indicator_reactions is not None and problem_type not in MILP_PROBLEM_TYPES:
        raise SpectraValueError(
            f"indicator_reactions needs a mixed-integer formulation, one of "
            f"{MILP_PROBLEM_TYPES}, not {problem_type!r}."
        )
    if indicator_groups is not None and problem_type not in MILP_PROBLEM_TYPES:
        raise SpectraValueError(
            f"indicator_groups needs a mixed-integer formulation, one of "
            f"{MILP_PROBLEM_TYPES}, not {problem_type!r}."
        )

    if problem_type == MIN_NET_LP:
        kept, flux_map, n_lps = _min_net_lp_core(
            model=model,
            directions=directions,
            weights=weights,
            tol=tol,
            steady_state=steady_state,
            inclusion_cutoff=inclusion_cutoff,
            feasibility_tolerance=feasibility_tolerance,
            optimality_tolerance=optimality_tolerance,
            run_irls=True,
            pass0_flux_ready=pass0_flux_ready,
        )
        v_arr = np.array([flux_map[r.id] for r in model.reactions], dtype=float)
        return MilpSolution(included=kept, selected=kept), v_arr, n_lps

    if problem_type == OPTIM_BIOMASS:
        kept, flux_map, n_lps = _optim_biomass_core(
            model=model,
            directions=directions,
            weights=weights,
            tol=tol,
            steady_state=steady_state,
            inclusion_cutoff=inclusion_cutoff,
            feasibility_tolerance=feasibility_tolerance,
            optimality_tolerance=optimality_tolerance,
        )
        v_arr = np.array([flux_map[r.id] for r in model.reactions], dtype=float)
        return MilpSolution(included=kept, selected=kept), v_arr, n_lps

    if problem_type == MIN_NET_MILP:
        sol, v_arr = _min_net_milp_core(
            model=model,
            directions=directions,
            weights=weights,
            tol=tol,
            steady_state=steady_state,
            time_limit=time_limit,
            previous_solutions=previous_solutions,
            indicator_reactions=indicator_reactions,
            inclusion_cutoff=inclusion_cutoff,
            indicator_groups=indicator_groups,
            feasibility_tolerance=feasibility_tolerance,
            optimality_tolerance=optimality_tolerance,
            integrality_tolerance=integrality_tolerance,
            warm_start=warm_start,
            auto_warm_start=True,
        )
        return sol, v_arr, 0

    sol, v_arr = _trade_off_core(
        model=model,
        directions=directions,
        weights=weights,
        tol=tol,
        steady_state=steady_state,
        time_limit=time_limit,
        previous_solutions=previous_solutions,
        indicator_reactions=indicator_reactions,
        inclusion_cutoff=inclusion_cutoff,
        indicator_groups=indicator_groups,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        integrality_tolerance=integrality_tolerance,
        warm_start=warm_start,
        auto_warm_start=True,
    )
    return sol, v_arr, 0


def _extract(model: "Model", keep_ids: Set[str], remove_genes: bool) -> "Model":
    """Return a copy of the model holding only the given reactions."""
    return subset_model(model, keep_ids, remove_genes)


def spectra_me(
    model: "Model",
    core_reactions: Iterable,
    tol: float = 1e-4,
    consistency_type: str = STOICHIOMETRY,
    weights: Optional[Dict[str, float]] = None,
    n_solutions: int = 1,
    alt_solution_method: Optional[str] = None,
    problem_type: str = MIN_NET_LP,
    time_limit: Optional[float] = 7200.0,
    remove_genes: bool = True,
    previous_solutions: Optional[List[Set[str]]] = None,
    seed: Optional[int] = None,
    inclusion_cutoff: Optional[float] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
    indicator_groups: Optional[Mapping[str, Iterable[str]]] = None,
    return_solutions: bool = False,
    feasibility_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    optimality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    integrality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
) -> Union["Model", List["Model"], Tuple]:
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
        The minimum absolute flux every core reaction in the extracted model
        has to carry (default 1e-4).
    consistency_type : {"stoichiometry", "topology"}, optional
        Whether to assume a steady state (``S v = 0``) or an accumulation
        condition (``S v >= 0``) (default "stoichiometry").
    weights : dict of {str: float}, optional
        The weight of each reaction (or group), keyed by identifier; reactions
        left out get 1.0. For ``"minNetLP"``, ``"minNetMILP"``, and
        ``"optimBiomass"`` weights must be non-negative costs; ``"tradeOff"``
        accepts any finite real number (positive rewards, zero neutral,
        negative penalties; default None, i.e. all ones).
    n_solutions : int, optional
        How many alternative models to return (default 1).
    alt_solution_method : {"coreDirection", "pathwayExclusion"}, optional
        How to find alternatives when ``n_solutions > 1``. When ``None`` (the
        default), selects ``"coreDirection"`` for LP formulations
        (``"minNetLP"``, ``"optimBiomass"``) and ``"pathwayExclusion"`` for
        MILP formulations (``"minNetMILP"``, ``"tradeOff"``).
    problem_type : {"minNetLP", "minNetMILP", "tradeOff", "optimBiomass"}, optional
        The network inference formulation, one of :data:`PROBLEM_TYPES`
        (default "minNetLP").
    time_limit : float, optional
        The maximum time to spend on each mixed-integer solve, in seconds
        (default 7200).
    remove_genes : bool, optional
        Whether to drop the genes left without a reaction by the extraction
        (default True). Orphaned metabolites are always dropped.
    previous_solutions : list of set of str, optional
        Indicator sets to exclude from the first solution too, for a
        mixed-integer `problem_type` (default None).
    seed : int, optional
        A seed for the random coefficients, making the result reproducible
        (default None).
    inclusion_cutoff : float, optional
        Minimum absolute flux to retain a non-core reaction (default derived
        from ``(tol, feasibility_tolerance)`` for LP formulations and
        ``tol * 1e-7`` for MILP formulations).
    indicator_reactions : iterable of str, optional
        Restrict the mixed-integer binaries to these reactions; the rest
        stay continuous and are kept whatever they do.
    indicator_groups : dict of {str: iterable of str}, optional
        Give the reactions of each group one **shared** binary, for
        ``"minNetMILP"`` or ``"tradeOff"``. Mutually exclusive with
        `indicator_reactions`.
    return_solutions : bool, optional
        Also return the :class:`~spectra_cobra.MilpSolution` behind each
        model (default False).
    feasibility_tolerance : float, optional
        Primal feasibility tolerance (default 1e-9). Pass None to leave the
        solver's own in place.
    optimality_tolerance : float, optional
        Dual feasibility tolerance (default 1e-9). Pass None to leave the
        solver's own in place.
    integrality_tolerance : float, optional
        Integrality tolerance for MILP binaries (default 1e-9). Pass None to
        leave the solver's own in place.

    Returns
    -------
    cobra.Model or list of cobra.Model
        The extracted model, or a list of models when ``n_solutions > 1``.
        With ``return_solutions=True``, a tuple of that and the matching
        :class:`~spectra_cobra.MilpSolution` (or list of them).

    Raises
    ------
    SpectraInfeasibleCoreError
        If the core reactions cannot all carry at least `tol` of flux at once.
    SpectraValueError
        If an argument value (`problem_type`, `weights`, `tol`, etc.) is invalid.

    See Also
    --------
    spectra_ccme : Extraction from a model that need not be flux consistent.
    spectra_me2 : Thermodynamic cone extraction without finite flux bounds.

    """
    _validate_positive_finite(tol, "tol")
    _warn_on_loose_tolerance(model, tol, feasibility_tolerance)
    models, solutions = _spectra_me(
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
        indicator_groups=indicator_groups,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        integrality_tolerance=integrality_tolerance,
    )
    extracted = models[0] if n_solutions == 1 else models
    if not return_solutions:
        return extracted
    return extracted, (solutions[0] if n_solutions == 1 else solutions)


def _spectra_me(
    model: "Model",
    core_reactions: Iterable,
    tol: float,
    consistency_type: str,
    weights: Optional[Dict[str, float]],
    n_solutions: int,
    alt_solution_method: Optional[str],
    problem_type: str,
    time_limit: Optional[float],
    remove_genes: bool,
    previous_solutions: Optional[List[Set[str]]],
    seed: Optional[int],
    blocked_ids: Optional[Set[str]] = None,
    inclusion_cutoff: Optional[float] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
    indicator_groups: Optional[Mapping[str, Iterable[str]]] = None,
    core_is_known_live: bool = False,
    feasibility_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    optimality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    integrality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
) -> Tuple[List["Model"], List[MilpSolution]]:
    """Run Stage-1 direction discovery and Stage-2 extraction for `spectra_me`."""
    if isinstance(n_solutions, bool) or n_solutions < 1:
        raise SpectraValueError(f"n_solutions must be at least 1, not {n_solutions!r}.")

    steady_state = validate_consistency_type(consistency_type)
    _validate_problem_type(problem_type)
    resolved_alt = _resolve_alt_solution_method(
        alt_solution_method, problem_type, n_solutions
    )

    if previous_solutions and problem_type not in MILP_PROBLEM_TYPES:
        raise SpectraValueError(
            f"Excluding previous solutions needs a mixed-integer formulation, "
            f"one of {MILP_PROBLEM_TYPES}, not {problem_type!r}."
        )
    if indicator_reactions is not None and problem_type not in MILP_PROBLEM_TYPES:
        raise SpectraValueError(
            f"indicator_reactions needs a mixed-integer formulation, one of "
            f"{MILP_PROBLEM_TYPES}, not {problem_type!r}."
        )
    if indicator_groups is not None and problem_type not in MILP_PROBLEM_TYPES:
        raise SpectraValueError(
            f"indicator_groups needs a mixed-integer formulation, one of "
            f"{MILP_PROBLEM_TYPES}, not {problem_type!r}."
        )

    tol_cfg = _resolve_lp_tolerances(
        model, tol, feasibility_tolerance, inclusion_cutoff
    )
    _apply_tolerances(
        model,
        feasibility_tolerance,
        optimality_tolerance,
        integrality_tolerance,
    )

    all_weights = _normalise_weights(
        model, weights, problem_type=problem_type, indicator_groups=indicator_groups
    )
    core_ids = _normalise_core(model, core_reactions)
    if blocked_ids:
        core_ids -= blocked_ids

    signs = reaction_signs(model)
    rng = np.random.default_rng(seed)
    has_per_rxn_weights = (indicator_groups is None) and (indicator_reactions is None)

    exclude = [set(s) for s in previous_solutions] if previous_solutions else []
    models: List["Model"] = []
    found: List[MilpSolution] = []

    # Stage 1: Primary core direction discovery (j = 1)
    directions, witness_flux, pass0_flux_ready, _ = _discover_core_directions(
        model=model,
        core_ids=core_ids,
        signs=signs,
        all_weights=all_weights,
        problem_type=problem_type,
        tol=tol_cfg.tol,
        steady_state=steady_state,
        inclusion_cutoff=inclusion_cutoff,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        lock_tol=tol_cfg.lock_tol,
        rng=rng,
        is_alt_sol=False,
        has_per_rxn_weights=has_per_rxn_weights,
        core_is_known_live=core_is_known_live,
    )

    try:
        solution, _, _ = _solve_formulation(
            model=model,
            directions=directions,
            weights=all_weights,
            tol=tol_cfg.tol,
            steady_state=steady_state,
            problem_type=problem_type,
            time_limit=time_limit,
            previous_solutions=exclude or None,
            inclusion_cutoff=inclusion_cutoff,
            indicator_reactions=indicator_reactions,
            indicator_groups=indicator_groups,
            feasibility_tolerance=feasibility_tolerance,
            optimality_tolerance=optimality_tolerance,
            integrality_tolerance=integrality_tolerance,
            warm_start=None,
            pass0_flux_ready=pass0_flux_ready,
        )
    except SpectraSolverError as err:
        if core_ids and not exclude:
            raise SpectraInfeasibleCoreError(
                f"No flux distribution lets all {len(core_ids)} core reactions "
                f"carry at least tol={tol_cfg.tol:.3g} at once. Use spectra_ccme "
                f"if the input model may contain blocked reactions."
            ) from err
        raise

    models.append(_extract(model, solution.included, remove_genes))
    found.append(solution)

    if n_solutions == 1:
        return models, found

    if resolved_alt == CORE_DIRECTION:
        can_both = any(
            oriented_bounds(model.reactions.get_by_id(rid), signs[rid])[0]
            <= -tol_cfg.lock_tol
            and oriented_bounds(model.reactions.get_by_id(rid), signs[rid])[1]
            >= tol_cfg.lock_tol
            for rid in core_ids
        )
        if not can_both:
            for _ in range(1, n_solutions):
                models.append(_extract(model, solution.included, remove_genes))
                found.append(solution)
            return models, found

        for _ in range(1, n_solutions):
            dir_j, _, _, _ = _discover_core_directions(
                model=model,
                core_ids=core_ids,
                signs=signs,
                all_weights=all_weights,
                problem_type=problem_type,
                tol=tol_cfg.tol,
                steady_state=steady_state,
                inclusion_cutoff=inclusion_cutoff,
                feasibility_tolerance=feasibility_tolerance,
                optimality_tolerance=optimality_tolerance,
                lock_tol=tol_cfg.lock_tol,
                rng=rng,
                is_alt_sol=True,
                has_per_rxn_weights=has_per_rxn_weights,
                core_is_known_live=core_is_known_live,
            )
            try:
                sol_j, _, _ = _solve_formulation(
                    model=model,
                    directions=dir_j,
                    weights=all_weights,
                    tol=tol_cfg.tol,
                    steady_state=steady_state,
                    problem_type=problem_type,
                    time_limit=time_limit,
                    previous_solutions=None,
                    inclusion_cutoff=inclusion_cutoff,
                    indicator_reactions=indicator_reactions,
                    indicator_groups=indicator_groups,
                    feasibility_tolerance=feasibility_tolerance,
                    optimality_tolerance=optimality_tolerance,
                    integrality_tolerance=integrality_tolerance,
                    warm_start=None,
                    pass0_flux_ready=None,
                )
            except SpectraSolverError as err:
                if core_ids:
                    raise SpectraInfeasibleCoreError(
                        f"No flux distribution lets all {len(core_ids)} core "
                        f"reactions carry at least tol={tol_cfg.tol:.3g} at once."
                    ) from err
                raise
            models.append(_extract(model, sol_j.included, remove_genes))
            found.append(sol_j)
        return models, found

    # pathwayExclusion
    for index in range(1, n_solutions):
        if not solution.selected:
            logger.info(
                "No non-empty indicator selection to exclude; returning %d "
                "solution(s) rather than %d.",
                index,
                n_solutions,
            )
            break
        exclude.append(solution.selected)
        try:
            solution, _, _ = _solve_formulation(
                model=model,
                directions=directions,
                weights=all_weights,
                tol=tol_cfg.tol,
                steady_state=steady_state,
                problem_type=problem_type,
                time_limit=time_limit,
                previous_solutions=exclude,
                inclusion_cutoff=inclusion_cutoff,
                indicator_reactions=indicator_reactions,
                indicator_groups=indicator_groups,
                feasibility_tolerance=feasibility_tolerance,
                optimality_tolerance=optimality_tolerance,
                integrality_tolerance=integrality_tolerance,
                warm_start=None,
                pass0_flux_ready=None,
            )
        except SpectraSolverError:
            logger.info(
                "No further alternative solution exists; returning the %d "
                "found rather than the %d requested.",
                index,
                n_solutions,
            )
            break
        models.append(_extract(model, solution.included, remove_genes))
        found.append(solution)

    return models, found


# =========================================================================
# spectra_ccme: Core-Scoped Consistency + Anchor-Refined Extraction
# =========================================================================


def _run_pinned_cc_loop(
    model: "Model",
    init_directions: Mapping[str, int],
    init_u: Iterable[str],
    signs: Mapping[str, float],
    can_fwd: Set[str],
    can_rev: Set[str],
    tol: float,
    lock_tol: float,
    rng: np.random.Generator,
    rand_order: bool = False,
) -> Tuple[Dict[str, int], Optional[Dict[str, float]], Set[str], int]:
    """Run core-scoped iterative `forward_cc`/`reverse` with already-settled core pinned."""
    directions = dict(init_directions)
    u_set = set(init_u)
    last_flux: Optional[Dict[str, float]] = None
    n_lps = 0

    prev_batch: Set[str] = {r for r, d in directions.items() if d != 0}
    prev_batch_flux: Optional[Dict[str, float]] = None

    while u_set:
        before = len(u_set)
        first_is_fwd = bool(rng.integers(2) == 0) if rand_order else True

        ordered_u = sorted(u_set)
        iter_weights = dict(
            zip(ordered_u, rng.uniform(*WEIGHT_RANGE, size=len(ordered_u)))
        )

        pass_order = (first_is_fwd, not first_is_fwd)
        for do_fwd in pass_order:
            if not u_set:
                break
            target = (u_set & can_fwd) if do_fwd else (u_set & can_rev)
            if not target:
                continue

            def _solve_cc_pass() -> Optional[Dict[str, float]]:
                with model:
                    apply_direction_bounds(model, directions, dict(signs), tol)
                    if do_fwd:
                        return forward_cc(
                            model,
                            target,
                            dict(signs),
                            tol,
                            rng,
                            objective_weights=iter_weights,
                        )
                    return reverse(
                        model,
                        target,
                        dict(signs),
                        tol,
                        rng,
                        objective_weights=iter_weights,
                    )

            f = _solve_cc_pass()
            n_lps += 1

            if f is None:
                if prev_batch_flux is not None and prev_batch:
                    borderline = {
                        r for r in prev_batch if abs(prev_batch_flux.get(r, 0.0)) < tol
                    }
                    if borderline:
                        for r in borderline:
                            directions[r] = 0
                        f = _solve_cc_pass()
                        n_lps += 1
                if f is None:
                    return directions, None, u_set, n_lps

            last_flux = f
            newly = {r for r in u_set if abs(f[r]) >= lock_tol}
            if newly:
                for r in newly:
                    ori = signs[r] * f[r]
                    directions[r] = 1 if ori > 0.0 else -1
                u_set -= newly
                prev_batch = newly
                prev_batch_flux = f

        if len(u_set) == before:
            break

    return directions, last_flux, u_set, n_lps


def _verify_blocked_core_unpinned(
    model: "Model",
    cand_u: Iterable[str],
    signs: Mapping[str, float],
    can_fwd: Set[str],
    can_rev: Set[str],
    tol: float,
    lock_tol: float,
    rng: np.random.Generator,
) -> Tuple[Set[str], Dict[str, int], int]:
    """Verify whether unsettled core reactions are individually blocked when unpinned."""
    u_set = set(cand_u)
    signs_found: Dict[str, int] = {}
    n_lps = 0

    while u_set:
        before = len(u_set)
        ordered_u = sorted(u_set)
        iter_weights = dict(
            zip(ordered_u, rng.uniform(*WEIGHT_RANGE, size=len(ordered_u)))
        )

        fwd_target = u_set & can_fwd
        if fwd_target:
            with model:
                f1 = forward_cc(
                    model,
                    fwd_target,
                    dict(signs),
                    tol,
                    rng,
                    objective_weights=iter_weights,
                )
            n_lps += 1
            if f1 is not None:
                hit1 = {r for r in u_set if abs(f1[r]) >= lock_tol}
                for r in hit1:
                    signs_found[r] = 1 if (signs[r] * f1[r]) > 0.0 else -1
                u_set -= hit1

        rev_target = u_set & can_rev
        if rev_target:
            with model:
                f2 = reverse(
                    model,
                    rev_target,
                    dict(signs),
                    tol,
                    rng,
                    objective_weights=iter_weights,
                )
            n_lps += 1
            if f2 is not None:
                hit2 = {r for r in u_set if abs(f2[r]) >= lock_tol}
                for r in hit2:
                    signs_found[r] = 1 if (signs[r] * f2[r]) > 0.0 else -1
                u_set -= hit2

        if len(u_set) == before:
            break

    return u_set, signs_found, n_lps


def _discover_consistent_core_directions(
    model: "Model",
    core_ids: Set[str],
    signs: Mapping[str, float],
    s_csc: "object",
    orig_lb: np.ndarray,
    orig_ub: np.ndarray,
    rxn_ids: List[str],
    has_coupling: bool,
    steady_state: bool,
    tol_cfg: _LpToleranceConfig,
    rng: np.random.Generator,
    is_alt_sol: bool = False,
) -> Tuple[Dict[str, int], Optional[Dict[str, float]], Set[str], np.ndarray, int]:
    """Stage 1a: 0-LP dead-end pruning + core-scoped `forward_cc`/`reverse` + blocked-core guard."""
    n = len(rxn_ids)
    directions = {rid: 0 for rid in rxn_ids}
    last_flux: Optional[Dict[str, float]] = None
    truly_blocked_core: Set[str] = set()
    n_lps = 0

    open_cols = (np.abs(orig_lb) >= tol_cfg.dead_end_zero_tol) | (
        np.abs(orig_ub) >= tol_cfg.dead_end_zero_tol
    )
    if not has_coupling:
        surv_cols = _strip_dead_ends(
            s_csc,
            orig_lb,
            orig_ub,
            open_cols,
            np.zeros(n, dtype=bool),
            steady_state,
            tol_cfg.dead_end_zero_tol,
        )
    else:
        surv_cols = open_cols

    # Skip 1: Empty core set -> 0 Stage-1 LPs
    if not core_ids:
        return directions, last_flux, truly_blocked_core, surv_cols, n_lps

    can_fwd_all: Set[str] = set()
    can_rev_all: Set[str] = set()
    surv_set = {rxn_ids[int(i)] for i in np.flatnonzero(surv_cols)}
    for i, rid in enumerate(rxn_ids):
        if not surv_cols[i]:
            continue
        lo, hi = oriented_bounds(model.reactions[i], signs[rid])
        if hi >= tol_cfg.lock_tol:
            can_fwd_all.add(rid)
        if lo <= -tol_cfg.lock_tol:
            can_rev_all.add(rid)

    can_reach = can_fwd_all | can_rev_all
    if not is_alt_sol:
        truly_blocked_core = {
            rid for rid in core_ids if (rid not in surv_set) or (rid not in can_reach)
        }

    cand_core = {
        rid for rid in core_ids if (rid in surv_set) and (rid in can_reach)
    }
    if not cand_core:
        return directions, last_flux, truly_blocked_core, surv_cols, n_lps

    dead_indices = np.flatnonzero(~surv_cols)
    with model, relaxed_mass_balance(model, steady_state):
        for idx in dead_indices:
            model.reactions[int(idx)].bounds = (0.0, 0.0)

        directions, last_flux, rem_u, loop_lps = _run_pinned_cc_loop(
            model=model,
            init_directions=directions,
            init_u=cand_core,
            signs=signs,
            can_fwd=can_fwd_all,
            can_rev=can_rev_all,
            tol=tol_cfg.tol,
            lock_tol=tol_cfg.lock_tol,
            rng=rng,
            rand_order=is_alt_sol,
        )
        n_lps += loop_lps

        if not rem_u or is_alt_sol:
            return directions, last_flux, truly_blocked_core, surv_cols, n_lps

        # Unpinned verification guard on remaining unsettled core reactions
        still_blocked_u, unpinned_signs, unp_lps = _verify_blocked_core_unpinned(
            model=model,
            cand_u=rem_u,
            signs=signs,
            can_fwd=can_fwd_all,
            can_rev=can_rev_all,
            tol=tol_cfg.tol,
            lock_tol=tol_cfg.lock_tol,
            rng=rng,
        )
        n_lps += unp_lps

        truly_blocked_core |= still_blocked_u
        conflict_or_greedy = rem_u - still_blocked_u

        if conflict_or_greedy:
            seed_dir = {rid: 0 for rid in rxn_ids}
            for rid in conflict_or_greedy:
                seed_dir[rid] = unpinned_signs[rid]
            rest_core = (cand_core - still_blocked_u) - conflict_or_greedy
            dir_rescue, flux_rescue, _, resc_lps = _run_pinned_cc_loop(
                model=model,
                init_directions=seed_dir,
                init_u=rest_core,
                signs=signs,
                can_fwd=can_fwd_all,
                can_rev=can_rev_all,
                tol=tol_cfg.tol,
                lock_tol=tol_cfg.lock_tol,
                rng=rng,
                rand_order=False,
            )
            n_lps += resc_lps

            n_active_cur = sum(1 for d in directions.values() if d != 0)
            n_active_rescue = sum(1 for d in dir_rescue.values() if d != 0)
            if flux_rescue is not None and n_active_rescue > n_active_cur:
                directions = dir_rescue
                last_flux = flux_rescue
            elif flux_rescue is None and n_active_cur == 0:
                dir_rescue2, flux_rescue2, _, resc_lps2 = _run_pinned_cc_loop(
                    model=model,
                    init_directions={rid: 0 for rid in rxn_ids},
                    init_u=conflict_or_greedy,
                    signs=signs,
                    can_fwd=can_fwd_all,
                    can_rev=can_rev_all,
                    tol=tol_cfg.tol,
                    lock_tol=tol_cfg.lock_tol,
                    rng=rng,
                    rand_order=False,
                )
                n_lps += resc_lps2
                if sum(1 for d in dir_rescue2.values() if d != 0) > n_active_cur:
                    directions = dir_rescue2
                    last_flux = flux_rescue2

            live_cand = cand_core - still_blocked_u
            dropped_mutual = [
                rid for rid in rxn_ids if rid in live_cand and directions[rid] == 0
            ]
            if dropped_mutual:
                active_cnt = sum(1 for d in directions.values() if d != 0)
                sample = ", ".join(dropped_mutual[:5])
                msg = (
                    f"{len(dropped_mutual)} core reaction(s) are individually "
                    f"consistent at tol={tol_cfg.tol:g} (so not returned in "
                    f"blocked_core), but could not carry >= tol simultaneously "
                    f"with the other {active_cnt} active core reactions under "
                    f"finite model bounds (e.g., {sample}). Keeping maximal "
                    f"simultaneously consistent core subset."
                )
                warnings.warn(msg, UserWarning, stacklevel=3)
                logger.warning(msg)

    return directions, last_flux, truly_blocked_core, surv_cols, n_lps


def _refine_rev_core_via_anchor_pinned(
    work_model: "Model",
    active_core_ids: Set[str],
    two_way_act: Set[str],
    prev_directions: Mapping[str, int],
    signs: Mapping[str, float],
    can_fwd_work: Set[str],
    can_rev_work: Set[str],
    w_anchor: Mapping[str, float],
    tol: float,
    steady_state: bool,
    inclusion_cutoff: Optional[float],
    feasibility_tolerance: Optional[float],
    optimality_tolerance: Optional[float],
    tol_cfg: _LpToleranceConfig,
    rng: np.random.Generator,
) -> Tuple[Dict[str, int], Optional[Dict[str, float]], int, bool, bool]:
    """Stage 1b: Refine reversible active core directions via 1-LP Anchor + `mPinned`."""
    exact_pass0 = False
    directions = dict(prev_directions)
    for rid in two_way_act:
        directions[rid] = 0

    try:
        _, v_anchor, n_lps = _min_net_lp_core(
            model=work_model,
            directions=directions,
            weights=w_anchor,
            tol=tol,
            steady_state=steady_state,
            inclusion_cutoff=inclusion_cutoff,
            feasibility_tolerance=feasibility_tolerance,
            optimality_tolerance=optimality_tolerance,
            run_irls=False,
        )
    except SpectraSolverError:
        return dict(prev_directions), None, 1, False, False

    last_flux: Optional[Dict[str, float]] = v_anchor
    for rid in two_way_act:
        if abs(v_anchor[rid]) >= DIRECTION_CUTOFF_FACTOR * tol:
            ori = signs[rid] * v_anchor[rid]
            directions[rid] = 1 if ori > 0.0 else -1

    unsettled_rev = {r for r in two_way_act if directions[r] == 0}
    if not unsettled_rev:
        if all(abs(v_anchor[r]) >= tol for r in active_core_ids):
            return directions, last_flux, n_lps, True, True
        if all(abs(v_anchor[r]) >= tol_cfg.lock_tol for r in active_core_ids):
            return directions, last_flux, n_lps, True, False
        for rid in two_way_act:
            if abs(v_anchor[rid]) < tol_cfg.lock_tol:
                directions[rid] = 0
        unsettled_rev = {r for r in two_way_act if directions[r] == 0}

    with work_model, relaxed_mass_balance(work_model, steady_state):
        directions, flux_pin, rem_u, pin_lps = _run_pinned_cc_loop(
            model=work_model,
            init_directions=directions,
            init_u=unsettled_rev,
            signs=signs,
            can_fwd=can_fwd_work,
            can_rev=can_rev_work,
            tol=tol,
            lock_tol=tol_cfg.lock_tol,
            rng=rng,
            rand_order=False,
        )
    n_lps += pin_lps
    if flux_pin is not None:
        last_flux = flux_pin
    ok_all = (not rem_u) and (
        sum(1 for d in directions.values() if d != 0) == len(active_core_ids)
    )
    return directions, last_flux, n_lps, ok_all, exact_pass0


def _run_ccme_single_extraction(
    work_model: "Model",
    directions: Dict[str, int],
    all_weights: Dict[str, float],
    problem_type: str,
    tol: float,
    steady_state: bool,
    time_limit: Optional[float],
    witness_flux: Optional[Mapping[str, float]],
    pass0_flux_ready: Optional[Mapping[str, float]],
    previous_solutions: Optional[List[Set[str]]],
    indicator_reactions: Optional[Iterable[str]],
    indicator_groups: Optional[Mapping[str, Iterable[str]]],
    inclusion_cutoff: Optional[float],
    feasibility_tolerance: Optional[float],
    optimality_tolerance: Optional[float],
    integrality_tolerance: Optional[float],
    tol_cfg: _LpToleranceConfig,
) -> MilpSolution:
    """Run Stage-2 extraction + consistency certification on `work_model` in `spectra_ccme`."""
    rxn_ids = [r.id for r in work_model.reactions]
    core_pinned_mask = np.array(
        [directions.get(rid, 0) != 0 for rid in rxn_ids], dtype=bool
    )

    if problem_type == MIN_NET_LP:
        try:
            sol, _, _ = _solve_formulation(
                model=work_model,
                directions=directions,
                weights=all_weights,
                tol=tol,
                steady_state=steady_state,
                problem_type=MIN_NET_LP,
                time_limit=time_limit,
                inclusion_cutoff=inclusion_cutoff,
                feasibility_tolerance=feasibility_tolerance,
                optimality_tolerance=optimality_tolerance,
                integrality_tolerance=integrality_tolerance,
                pass0_flux_ready=pass0_flux_ready,
            )
            return sol
        except SpectraSolverError:
            if witness_flux is not None:
                s_csc, w_lb, w_ub, _, has_cpl, extra_c = _build_network_arrays(
                    work_model
                )
                v_wit = np.array(
                    [float(witness_flux.get(rid, 0.0)) for rid in rxn_ids], dtype=float
                )
                raw_inc = {
                    rxn_ids[i]
                    for i in range(len(rxn_ids))
                    if core_pinned_mask[i] or abs(v_wit[i]) >= tol_cfg.supp_tol
                }
                cert_inc, _ = _certify_milp_support(
                    model=work_model,
                    s_csc=s_csc,
                    lb=w_lb,
                    ub=w_ub,
                    rxn_ids=rxn_ids,
                    raw_included=raw_inc,
                    v_sol=v_wit,
                    core_pinned_mask=core_pinned_mask,
                    has_groups=False,
                    steady_state=steady_state,
                    has_coupling=has_cpl,
                    extra_cons=extra_c,
                    feasibility_tolerance=feasibility_tolerance,
                    optimality_tolerance=optimality_tolerance,
                    tol_cfg=tol_cfg,
                )
                return MilpSolution(included=cert_inc, selected=cert_inc)
            raise

    if problem_type == OPTIM_BIOMASS:
        sol, _, _ = _solve_formulation(
            model=work_model,
            directions=directions,
            weights=all_weights,
            tol=tol,
            steady_state=steady_state,
            problem_type=OPTIM_BIOMASS,
            time_limit=time_limit,
            inclusion_cutoff=inclusion_cutoff,
            feasibility_tolerance=feasibility_tolerance,
            optimality_tolerance=optimality_tolerance,
            integrality_tolerance=integrality_tolerance,
        )
        return sol

    # MILP formulations: minNetMILP or tradeOff
    sol, v_arr, _ = _solve_formulation(
        model=work_model,
        directions=directions,
        weights=all_weights,
        tol=tol,
        steady_state=steady_state,
        problem_type=problem_type,
        time_limit=time_limit,
        previous_solutions=previous_solutions,
        inclusion_cutoff=inclusion_cutoff,
        indicator_reactions=indicator_reactions,
        indicator_groups=indicator_groups,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        integrality_tolerance=integrality_tolerance,
        warm_start=None,
    )
    if not sol.included:
        return sol

    s_csc, w_lb, w_ub, _, has_cpl, extra_c = _build_network_arrays(work_model)
    has_groups = (indicator_groups is not None) or (indicator_reactions is not None)
    cert_inc, _ = _certify_milp_support(
        model=work_model,
        s_csc=s_csc,
        lb=w_lb,
        ub=w_ub,
        rxn_ids=rxn_ids,
        raw_included=sol.included,
        v_sol=v_arr,
        core_pinned_mask=core_pinned_mask,
        has_groups=has_groups,
        steady_state=steady_state,
        has_coupling=has_cpl,
        extra_cons=extra_c,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        tol_cfg=tol_cfg,
    )
    return MilpSolution(included=cert_inc, selected=sol.selected)


def spectra_ccme(
    model: "Model",
    core_reactions: Iterable,
    tol: float = 1e-4,
    consistency_type: str = STOICHIOMETRY,
    weights: Optional[Dict[str, float]] = None,
    n_solutions: int = 1,
    alt_solution_method: Optional[str] = None,
    problem_type: str = MIN_NET_LP,
    time_limit: Optional[float] = 7200.0,
    remove_genes: bool = True,
    seed: Optional[int] = None,
    inclusion_cutoff: Optional[float] = None,
    feasibility_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    optimality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    integrality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    indicator_reactions: Optional[Iterable[str]] = None,
    indicator_groups: Optional[Mapping[str, Iterable[str]]] = None,
    return_solutions: bool = False,
) -> Tuple:
    """Check consistency and extract a model in a single pass.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from. It need not be flux consistent.
    core_reactions : iterable
        The reactions requested in the extracted model, as identifiers,
        :class:`cobra.Reaction` objects, or indices. Blocked core reactions
        are omitted and returned in `blocked_core` rather than raising.
    tol : float, optional
        The minimum absolute flux required for consistency and core activity
        (default 1e-4).
    consistency_type : {"stoichiometry", "topology"}, optional
        Whether to assume a steady state (``S v = 0``) or an accumulation
        condition (``S v >= 0``) (default "stoichiometry").
    weights : dict of {str: float}, optional
        Reaction (or group) weights (default None, i.e. all ones). Must be
        non-negative for ``"minNetLP"``, ``"minNetMILP"``, and
        ``"optimBiomass"``; any finite real numbers for ``"tradeOff"``.
    n_solutions : int, optional
        How many alternative models to return (default 1).
    alt_solution_method : {"coreDirection", "pathwayExclusion"}, optional
        How to find alternatives when ``n_solutions > 1``. Defaults to
        ``"coreDirection"`` for LP formulations and ``"pathwayExclusion"``
        for MILP formulations.
    problem_type : {"minNetLP", "minNetMILP", "tradeOff", "optimBiomass"}, optional
        The network inference formulation, one of :data:`PROBLEM_TYPES`
        (default "minNetLP").
    time_limit : float, optional
        The maximum time to spend on each mixed-integer solve, in seconds
        (default 7200).
    remove_genes : bool, optional
        Whether to drop the genes left without a reaction by the extraction
        (default True). Orphaned metabolites are always dropped.
    seed : int, optional
        A seed for the random coefficients (default None).
    inclusion_cutoff : float, optional
        Minimum absolute flux to retain a non-core reaction (default derived
        from ``(tol, feasibility_tolerance)`` for LP formulations and
        ``tol * 1e-7`` for MILP formulations).
    feasibility_tolerance : float, optional
        Primal feasibility tolerance (default 1e-9).
    optimality_tolerance : float, optional
        Dual feasibility tolerance (default 1e-9).
    integrality_tolerance : float, optional
        Integrality tolerance for MILP binaries (default 1e-9).
    indicator_reactions : iterable of str, optional
        Restrict MILP binaries to these reactions.
    indicator_groups : dict of {str: iterable of str}, optional
        Shared MILP binary groups for ``"minNetMILP"`` or ``"tradeOff"``.
    return_solutions : bool, optional
        If True, append the :class:`~spectra_cobra.MilpSolution` (or list of
        them) as a third element of the returned tuple (default False).

    Returns
    -------
    tuple of (cobra.Model or list of cobra.Model, list of str)
        The extracted consistent model (or list of models when
        ``n_solutions > 1``) and the identifiers of the core reactions that
        are blocked at `tol` in `model`.

    See Also
    --------
    spectra_cc : The consistency check on its own.
    spectra_me : Extraction from a model already known to be consistent.

    """
    if isinstance(n_solutions, bool) or n_solutions < 1:
        raise SpectraValueError(f"n_solutions must be at least 1, not {n_solutions!r}.")

    steady_state = validate_consistency_type(consistency_type)
    _validate_problem_type(problem_type)
    resolved_alt = _resolve_alt_solution_method(
        alt_solution_method, problem_type, n_solutions
    )
    if indicator_reactions is not None and problem_type not in MILP_PROBLEM_TYPES:
        raise SpectraValueError(
            f"indicator_reactions needs a mixed-integer formulation, one of "
            f"{MILP_PROBLEM_TYPES}, not {problem_type!r}."
        )
    if indicator_groups is not None and problem_type not in MILP_PROBLEM_TYPES:
        raise SpectraValueError(
            f"indicator_groups needs a mixed-integer formulation, one of "
            f"{MILP_PROBLEM_TYPES}, not {problem_type!r}."
        )

    tol_cfg = _resolve_lp_tolerances(
        model, tol, feasibility_tolerance, inclusion_cutoff
    )
    _warn_on_loose_tolerance(model, tol_cfg.tol, feasibility_tolerance)
    _apply_tolerances(
        model,
        feasibility_tolerance,
        optimality_tolerance,
        integrality_tolerance,
    )

    all_weights = _normalise_weights(
        model, weights, problem_type=problem_type, indicator_groups=indicator_groups
    )
    core_ids = _normalise_core(model, core_reactions)
    signs = reaction_signs(model)
    rng = np.random.default_rng(seed)
    has_per_rxn_weights = (indicator_groups is None) and (indicator_reactions is None)

    s_csc, orig_lb, orig_ub, rxn_ids, has_coupling, _ = _build_network_arrays(model)

    # Stage 1a: 0-LP dead-end pruning + core-scoped forward_cc/reverse + unpinned guard
    (
        directions,
        witness_flux,
        truly_blocked_core,
        surv_cols,
        _,
    ) = _discover_consistent_core_directions(
        model=model,
        core_ids=core_ids,
        signs=signs,
        s_csc=s_csc,
        orig_lb=orig_lb,
        orig_ub=orig_ub,
        rxn_ids=rxn_ids,
        has_coupling=has_coupling,
        steady_state=steady_state,
        tol_cfg=tol_cfg,
        rng=rng,
        is_alt_sol=False,
    )

    blocked_core = [rid for rid in rxn_ids if rid in truly_blocked_core]
    if blocked_core:
        logger.warning(
            "%d of the %d core reactions are blocked and were dropped: %s",
            len(blocked_core),
            len(core_ids),
            blocked_core[:5],
        )

    # Return empty model immediately in 0 Stage-2 solves if 100% of non-empty core is blocked
    if core_ids and not any(d != 0 for d in directions.values()):
        empty_sol = MilpSolution(included=set(), selected=set())
        if n_solutions == 1:
            empty_mod = _extract(model, set(), remove_genes)
            if return_solutions:
                return empty_mod, blocked_core, empty_sol
            return empty_mod, blocked_core
        empty_mods = [_extract(model, set(), remove_genes) for _ in range(n_solutions)]
        if return_solutions:
            return empty_mods, blocked_core, [empty_sol for _ in range(n_solutions)]
        return empty_mods, blocked_core

    shut_indices = [
        i
        for i, rid in enumerate(rxn_ids)
        if (not surv_cols[i]) or (rid in truly_blocked_core)
    ]

    with model:
        for idx in shut_indices:
            model.reactions[idx].bounds = (0.0, 0.0)

        # Stage 1b: Anchor-Refined Reversible Active Core Sign Selection
        active_core_primary = {r for r, d in directions.items() if d != 0}
        pass0_flux_ready: Optional[Dict[str, float]] = None

        can_fwd_work: Set[str] = set()
        can_rev_work: Set[str] = set()
        for rxn in model.reactions:
            lo, hi = oriented_bounds(rxn, signs[rxn.id])
            if hi >= tol_cfg.lock_tol:
                can_fwd_work.add(rxn.id)
            if lo <= -tol_cfg.lock_tol:
                can_rev_work.add(rxn.id)

        two_way_act = {
            r for r in active_core_primary if (r in can_fwd_work) and (r in can_rev_work)
        }
        one_way_act = active_core_primary - two_way_act

        use_anchor_ref = False
        w_anchor_ref: Optional[Dict[str, float]] = None
        if (
            has_per_rxn_weights
            and two_way_act
            and (one_way_act or not _zero_flux_is_feasible(model, steady_state))
        ):
            if problem_type in (MIN_NET_LP, MIN_NET_MILP):
                w_anchor_ref = {
                    rid: (0.0 if rid in active_core_primary else float(all_weights[rid]))
                    for rid in rxn_ids
                }
                use_anchor_ref = True
            elif problem_type == TRADE_OFF and any(
                float(all_weights[rid]) < 0.0 for rid in rxn_ids
            ):
                w_anchor_ref = {}
                for rid in rxn_ids:
                    if rid in active_core_primary:
                        w_anchor_ref[rid] = 0.0
                    else:
                        wv = float(all_weights[rid])
                        w_anchor_ref[rid] = -wv if wv < 0.0 else 0.0
                use_anchor_ref = True

        if use_anchor_ref and w_anchor_ref is not None:
            (
                dir_ref,
                wit_ref,
                _,
                ok_ref,
                exact_pass0,
            ) = _refine_rev_core_via_anchor_pinned(
                work_model=model,
                active_core_ids=active_core_primary,
                two_way_act=two_way_act,
                prev_directions=directions,
                signs=signs,
                can_fwd_work=can_fwd_work,
                can_rev_work=can_rev_work,
                w_anchor=w_anchor_ref,
                tol=tol_cfg.tol,
                steady_state=steady_state,
                inclusion_cutoff=inclusion_cutoff,
                feasibility_tolerance=feasibility_tolerance,
                optimality_tolerance=optimality_tolerance,
                tol_cfg=tol_cfg,
                rng=rng,
            )
            if ok_ref:
                directions = dir_ref
                witness_flux = wit_ref
                if exact_pass0 and problem_type == MIN_NET_LP:
                    pass0_flux_ready = wit_ref

        # Stage 2: Primary Network Extraction (j = 1)
        solution = _run_ccme_single_extraction(
            work_model=model,
            directions=directions,
            all_weights=all_weights,
            problem_type=problem_type,
            tol=tol_cfg.tol,
            steady_state=steady_state,
            time_limit=time_limit,
            witness_flux=witness_flux,
            pass0_flux_ready=pass0_flux_ready,
            previous_solutions=None,
            indicator_reactions=indicator_reactions,
            indicator_groups=indicator_groups,
            inclusion_cutoff=inclusion_cutoff,
            feasibility_tolerance=feasibility_tolerance,
            optimality_tolerance=optimality_tolerance,
            integrality_tolerance=integrality_tolerance,
            tol_cfg=tol_cfg,
        )

        kept_sets: List[Set[str]] = [solution.included]
        solutions_list: List[MilpSolution] = [solution]

        if n_solutions > 1:
            active_core_mask = core_ids - truly_blocked_core
            if resolved_alt == CORE_DIRECTION:
                can_both_act = any(
                    (rid in can_fwd_work) and (rid in can_rev_work)
                    for rid in active_core_mask
                )
                if not can_both_act:
                    for _ in range(1, n_solutions):
                        kept_sets.append(solution.included)
                        solutions_list.append(solution)
                else:
                    w_s_csc, w_lb, w_ub, _, w_cpl, _ = _build_network_arrays(model)
                    for _ in range(1, n_solutions):
                        dir_j, wit_j, _, _, _ = _discover_consistent_core_directions(
                            model=model,
                            core_ids=active_core_mask,
                            signs=signs,
                            s_csc=w_s_csc,
                            orig_lb=w_lb,
                            orig_ub=w_ub,
                            rxn_ids=rxn_ids,
                            has_coupling=w_cpl,
                            steady_state=steady_state,
                            tol_cfg=tol_cfg,
                            rng=rng,
                            is_alt_sol=True,
                        )
                        sol_j = _run_ccme_single_extraction(
                            work_model=model,
                            directions=dir_j,
                            all_weights=all_weights,
                            problem_type=problem_type,
                            tol=tol_cfg.tol,
                            steady_state=steady_state,
                            time_limit=time_limit,
                            witness_flux=wit_j,
                            pass0_flux_ready=None,
                            previous_solutions=None,
                            indicator_reactions=indicator_reactions,
                            indicator_groups=indicator_groups,
                            inclusion_cutoff=inclusion_cutoff,
                            feasibility_tolerance=feasibility_tolerance,
                            optimality_tolerance=optimality_tolerance,
                            integrality_tolerance=integrality_tolerance,
                            tol_cfg=tol_cfg,
                        )
                        kept_sets.append(sol_j.included)
                        solutions_list.append(sol_j)
            else:
                exclude: List[Set[str]] = []
                cur_sol = solution
                for index in range(1, n_solutions):
                    if not cur_sol.selected:
                        break
                    exclude.append(cur_sol.selected)
                    try:
                        cur_sol = _run_ccme_single_extraction(
                            work_model=model,
                            directions=directions,
                            all_weights=all_weights,
                            problem_type=problem_type,
                            tol=tol_cfg.tol,
                            steady_state=steady_state,
                            time_limit=time_limit,
                            witness_flux=None,
                            pass0_flux_ready=None,
                            previous_solutions=exclude,
                            indicator_reactions=indicator_reactions,
                            indicator_groups=indicator_groups,
                            inclusion_cutoff=inclusion_cutoff,
                            feasibility_tolerance=feasibility_tolerance,
                            optimality_tolerance=optimality_tolerance,
                            integrality_tolerance=integrality_tolerance,
                            tol_cfg=tol_cfg,
                        )
                    except SpectraSolverError:
                        logger.info(
                            "No further alternative solution exists; returning "
                            "the %d found rather than the %d requested.",
                            index,
                            n_solutions,
                        )
                        break
                    kept_sets.append(cur_sol.included)
                    solutions_list.append(cur_sol)

    # Subset from original `model` (outside `with model:`) so original bounds and S are preserved
    models_out = [_extract(model, ks, remove_genes) for ks in kept_sets]
    if n_solutions == 1:
        if return_solutions:
            return models_out[0], blocked_core, solutions_list[0]
        return models_out[0], blocked_core
    if return_solutions:
        return models_out, blocked_core, solutions_list
    return models_out, blocked_core


# =========================================================================
# spectra_me2: Thermodynamic Cone Extraction (Unbounded Flux Magnitudes)
# =========================================================================


def spectra_me2(
    model: "Model",
    core_reactions: Iterable,
    tol: float = 1e-4,
    weights: Optional[Dict[str, float]] = None,
    n_solutions: int = 1,
    remove_genes: bool = True,
    seed: Optional[int] = None,
    feasibility_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    optimality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
) -> Union["Model", List["Model"]]:
    """Extract a thermodynamic subnetwork without finite box bounds (`spectraME2`).

    Replaces finite box bounds with the pure thermodynamic cone —
    ``(-inf, inf)`` on reversible reactions, ``(0, inf)`` on forward-only
    reactions, and ``(-inf, 0)`` on reverse-only reactions (``ub <= 0``) —
    then extracts around `core_reactions` via Stage-1 direction discovery and
    Stage-2 Split-Variable ``min_net_lp`` while preserving each reaction's
    original cobrapy orientation.

    Parameters
    ----------
    model : cobra.Model
        Flux-consistent model to extract from.
    core_reactions : iterable
        Core reactions required in the extracted model (identifiers,
        :class:`cobra.Reaction` objects, or indices).
    tol : float, optional
        Minimum absolute flux required for core reactions (default 1e-4).
    weights : dict of {str: float}, optional
        Non-negative cost weights for reactions (default None, i.e. all ones).
        Negative weights raise :class:`ValueError`.
    n_solutions : int, optional
        Number of alternative solutions via ``"coreDirection"`` (default 1).
    remove_genes : bool, optional
        Whether to drop unused genes from the extracted model (default True).
    seed : int, optional
        RNG seed for reproducibility (default None).
    feasibility_tolerance : float, optional
        LP primal feasibility tolerance (default 1e-9).
    optimality_tolerance : float, optional
        LP dual optimality tolerance (default 1e-9).

    Returns
    -------
    cobra.Model or list of cobra.Model
        The extracted thermodynamic model, or a list of `n_solutions` models.

    Raises
    ------
    SpectraInfeasibleCoreError
        If the core reactions cannot all carry at least `tol` of flux at once.
    SpectraValueError
        If `weights` contains negative or non-finite values, or `tol <= 0`.

    """
    if isinstance(n_solutions, bool) or n_solutions < 1:
        raise SpectraValueError(f"n_solutions must be at least 1, not {n_solutions!r}.")

    tol_cfg = _resolve_lp_tolerances(model, tol, feasibility_tolerance, None)
    _warn_on_loose_tolerance(model, tol_cfg.tol, feasibility_tolerance)

    all_weights = _normalise_weights(model, weights, problem_type=MIN_NET_LP)
    core_ids = _normalise_core(model, core_reactions)

    thermo_model = subset_model(
        model,
        {rxn.id for rxn in model.reactions},
        remove_genes=False,
        keep_orphan_metabolites=True,
    )
    _apply_tolerances(thermo_model, feasibility_tolerance, optimality_tolerance)

    rev_ids: Set[str] = set()
    for rxn in thermo_model.reactions:
        orig_lb = float(rxn.lower_bound)
        orig_ub = float(rxn.upper_bound)
        if orig_lb < 0.0 and orig_ub > 0.0:
            rxn.bounds = (-float("inf"), float("inf"))
            rev_ids.add(rxn.id)
        elif orig_ub <= 0.0 and orig_lb < 0.0:
            rxn.bounds = (-float("inf"), 0.0)
        elif orig_lb >= 0.0 and orig_ub > 0.0:
            rxn.bounds = (0.0, float("inf"))
        else:
            rxn.bounds = (0.0, 0.0)

    signs = reaction_signs(thermo_model)
    rng = np.random.default_rng(seed)

    directions, _, pass0_ready, _ = _discover_core_directions(
        model=thermo_model,
        core_ids=core_ids,
        signs=signs,
        all_weights=all_weights,
        problem_type=MIN_NET_LP,
        tol=tol_cfg.tol,
        steady_state=True,
        inclusion_cutoff=None,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        lock_tol=tol_cfg.lock_tol,
        rng=rng,
        is_alt_sol=False,
        has_per_rxn_weights=True,
    )
    kept, _, _ = _min_net_lp_core(
        model=thermo_model,
        directions=directions,
        weights=all_weights,
        tol=tol_cfg.tol,
        steady_state=True,
        inclusion_cutoff=None,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        run_irls=True,
        pass0_flux_ready=pass0_ready,
    )
    first = _extract(thermo_model, kept, remove_genes)
    if n_solutions == 1:
        return first

    models_out = [first]
    if not (core_ids & rev_ids):
        for _ in range(1, n_solutions):
            models_out.append(_extract(thermo_model, kept, remove_genes))
        return models_out

    for _ in range(1, n_solutions):
        dir_j, _, _, _ = _discover_core_directions(
            model=thermo_model,
            core_ids=core_ids,
            signs=signs,
            all_weights=all_weights,
            problem_type=MIN_NET_LP,
            tol=tol_cfg.tol,
            steady_state=True,
            inclusion_cutoff=None,
            feasibility_tolerance=feasibility_tolerance,
            optimality_tolerance=optimality_tolerance,
            lock_tol=tol_cfg.lock_tol,
            rng=rng,
            is_alt_sol=True,
            has_per_rxn_weights=True,
        )
        kept_j, _, _ = _min_net_lp_core(
            model=thermo_model,
            directions=dir_j,
            weights=all_weights,
            tol=tol_cfg.tol,
            steady_state=True,
            inclusion_cutoff=None,
            feasibility_tolerance=feasibility_tolerance,
            optimality_tolerance=optimality_tolerance,
            run_irls=True,
            pass0_flux_ready=None,
        )
        models_out.append(_extract(thermo_model, kept_j, remove_genes))
    return models_out
