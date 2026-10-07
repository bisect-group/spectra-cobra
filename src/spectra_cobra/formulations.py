"""The network inference formulations SPECTRA extracts models with.

Each function here takes a model, a per-reaction direction read off the
iterative phase, and per-reaction weights, and returns the reactions that
belong in the extracted model. They differ in what they optimise:

==================  ==========================================  ==============
Formulation         Objective                                   Weights
==================  ==========================================  ==============
``min_net_lp``      minimise weighted total absolute flux       non-negative
``min_net_milp``    minimise weighted reaction count            non-negative
``trade_off``       maximise weighted included reaction count   any real
``growth_optim``    maximise biomass minus weighted flux        non-negative
==================  ==========================================  ==============

A reaction's *direction* is 1 if it must carry positive oriented flux, -1 if
it must carry negative oriented flux, and 0 if it is free to be dropped. The
directed reactions are the core set; the rest are what the formulations choose
between.
"""

import math
from dataclasses import dataclass
from logging import getLogger
from typing import (
    TYPE_CHECKING,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
)

from cobra.util.solver import linear_reaction_coefficients
from optlang.interface import FEASIBLE, OPTIMAL, TIME_LIMIT
from optlang.symbolics import Zero, add

from ._orientation import (
    apply_direction_bounds,
    is_reversible,
    oriented_bounds,
    reaction_signs,
    relaxed_mass_balance,
)
from .exceptions import SpectraError, SpectraSolverError

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)

#: Statuses a MILP may stop at and still carry a usable solution: optimal,
#: or a feasible solution found before the time limit. ``TIME_LIMIT`` belongs
#: here because that is what Gurobi reports when it stops at `time_limit`
#: holding an incumbent, which is the whole point of setting a limit; it can
#: also be reported with no solution at all, so :func:`_has_solution` checks.
ACCEPTABLE_MILP_STATUSES = (OPTIMAL, FEASIBLE, TIME_LIMIT)

#: The fraction of `tol` below which a flux counts as zero when reading the
#: extracted reaction set off a solution.
INCLUSION_CUTOFF_FACTOR = 1e-7


def _check_status(model: "Model", description: str, milp: bool = False) -> str:
    """Raise unless the solver reached a status carrying a usable solution.

    Parameters
    ----------
    model : cobra.Model
        The model that was just optimized.
    description : str
        What was being solved, for the error message.
    milp : bool, optional
        Whether a merely feasible solution is acceptable, as it is for a MILP
        stopped at its time limit (default False).

    Returns
    -------
    str
        The solver status.

    Raises
    ------
    SpectraSolverError
        If the status carries no usable solution.

    """
    status = model.solver.status
    acceptable = ACCEPTABLE_MILP_STATUSES if milp else (OPTIMAL,)
    if status not in acceptable:
        raise SpectraSolverError(
            f"The {description} terminated with status {status!r} instead of "
            f"reaching {'an acceptable status' if milp else 'optimality'}.",
            status=status,
        )
    if status != OPTIMAL:
        if not _has_solution(model):
            raise SpectraSolverError(
                f"The {description} terminated with status {status!r} before "
                f"finding any feasible solution. Allow it more time, or use a "
                f"linear formulation such as minNetLP.",
                status=status,
            )
        logger.warning(
            "The %s stopped at status %r rather than optimality; the returned "
            "model is feasible but may not be minimal.",
            description,
            status,
        )
    return status


def _has_solution(model: "Model") -> bool:
    """Return whether the solver is holding a solution that can be read.

    Parameters
    ----------
    model : cobra.Model
        The model that was just optimized.

    Returns
    -------
    bool
        Whether an objective value, and so a solution, is available.

    Notes
    -----
    A mixed-integer solve stopped at its time limit may or may not have found
    an incumbent, and the distinction is only visible by trying to read one.
    Each interface raises its own exception type when there is nothing to
    read, none of which can be imported without depending on that solver, so
    the exception is caught broadly on purpose.

    """
    try:
        value = model.solver.objective.value
    except Exception:  # noqa: BLE001 - see the note above.
        return False
    return value is not None and not math.isnan(value)


def _included_reactions(
    model: "Model", tol: float, inclusion_cutoff: Optional[float] = None
) -> Set[str]:
    """Return the reactions carrying flux in the current solution.

    Parameters
    ----------
    model : cobra.Model
        A model whose solver holds a solution.
    tol : float
        The flux threshold the formulation was given.
    inclusion_cutoff : float, optional
        The absolute flux at which a reaction counts as part of the model
        (default ``tol * 1e-7``).

    Returns
    -------
    set of str
        The identifiers of the reactions to keep.

    Notes
    -----
    The default is deliberately far below the solver's own tolerance, and
    that is not an oversight. Raising it to something "sensible" such as
    ``model.tolerance`` discards reactions that carry almost no flux but are
    nonetheless load-bearing for a mass balance, and the chains that depended
    on them then collapse: on Recon3D that leaves core reactions present in
    the extracted model but unable to carry flux. Keeping the cutoff low
    keeps those reactions, at the cost of also keeping some that are there
    only because of numerical noise. The trade is deliberate — a slightly
    dirtier model in exchange for core reactions that actually work.

    """
    if inclusion_cutoff is None:
        inclusion_cutoff = tol * INCLUSION_CUTOFF_FACTOR
    included = {rxn.id for rxn in model.reactions if abs(rxn.flux) >= inclusion_cutoff}
    _warn_on_noisy_support(model, included, tol)
    return included


def _warn_on_noisy_support(model: "Model", included: Set[str], tol: float) -> None:
    """Warn when the solution's support is partly below solver precision.

    Parameters
    ----------
    model : cobra.Model
        The solved model.
    included : set of str
        The reactions the solution is taken to include.
    tol : float
        The flux threshold the formulation was given.

    Notes
    -----
    A reaction carrying less flux than the solver's own feasibility
    tolerance is carrying a number the solver cannot distinguish from
    zero, so its presence in the support says nothing. When enough of
    the support is like that, the support is not a network: restricting
    the model to it leaves everything blocked, core reactions included.

    This separates cleanly rather than gradually, which is why it is
    worth reporting at all. On a two-tissue body of 18,090 reactions at
    ``tol=1e-4`` with ``model.tolerance=1e-9``, 645 of the 1,269 included
    reactions were below tolerance and 1,267 of them could not carry
    flux. At ``tol=1e-3`` the count was zero and the result verified.
    The cause is the gap between `tol` and the solver's tolerance being
    too small for the size of the problem, so the remedy is a larger
    `tol` -- not a different cutoff, which cannot recover information the
    solution never carried.

    """
    precision = model.tolerance
    if not precision:
        return
    noisy = sum(
        1 for rxn in model.reactions if rxn.id in included and abs(rxn.flux) < precision
    )
    if not noisy:
        return
    logger.warning(
        "%d of %d reactions in this solution carry less flux than the "
        "solver can resolve (model.tolerance=%.3g), so whether they belong "
        "in the model is not something the solution settles. The extracted "
        "model may be unable to carry flux at all, core reactions included. "
        "Raise tol (it is %.3g; ten times that usually suffices at this "
        "size), or use minNetMILP, whose membership comes from its binaries "
        "rather than from flux. check_extraction will tell you either way.",
        noisy,
        len(included),
        precision,
        tol,
    )


def _included_from_indicators(
    directions: Dict[str, int],
    indicators: Dict[str, "object"],
    groups: Mapping[str, Sequence[str]],
) -> Set[str]:
    """Return the reactions a mixed-integer solution keeps.

    Parameters
    ----------
    directions : dict of {str: int}
        The oriented direction of each reaction. The directed ones are kept
        unconditionally, since they are forced to carry flux.
    indicators : dict of {str: optlang.Variable}
        The binary inclusion variable of each indicator group.
    groups : dict of {str: sequence of str}
        The reactions each binary governs. One binary may govern several,
        in which case switching it on keeps all of them.

    Returns
    -------
    set of str
        The identifiers of the reactions to keep.

    Notes
    -----
    Reading the answer off the binaries rather than off the flux is both
    exact and necessary. A binary that the solver leaves at its integrality
    tolerance rather than exactly zero still multiplies a flux bound of up to
    a thousand, which leaks far more apparent flux than any sensible cutoff,
    so the flux-derived set can disagree with the binaries the objective and
    the exclusion constraints are actually written over.

    """
    kept = {rxn_id for rxn_id, direction in directions.items() if direction != 0}
    for key, var in indicators.items():
        if var.primal > 0.5:
            kept.update(groups[key])
    return kept


@dataclass(frozen=True)
class MilpSolution:
    """What a mixed-integer formulation found.

    Parameters
    ----------
    included : set of str
        The reactions the extracted model should contain.
    selected : set of str
        The **indicator groups** the solver switched on, named by their
        keys. With one binary per reaction -- the usual case -- those keys
        are reaction identifiers and this reads as the reactions selected.
        This is what :func:`_add_exclusion_constraints` has to be written
        over, and it is not always the same as `included`.

    Notes
    -----
    The two differ because a binary at zero does not pin its flux to zero
    exactly: :math:`v_i \\le ub_i z_i` with :math:`ub_i = 1000` still permits
    :math:`|v_i| \\le 10^{-6}` at the solver's integrality tolerance. A
    reaction carrying a trace flux the solution genuinely depends on can
    therefore be reported as switched off, and dropping it breaks the chain
    it was part of. `included` keeps anything carrying real flux whatever
    its indicator says; `selected` reports the indicators as they stand, so
    that an exclusion constraint forbids the assignment that actually
    occurred.

    """

    included: Set[str]
    selected: Set[str]


def _milp_solution(
    model: "Model",
    directions: Dict[str, int],
    indicators: Dict[str, "object"],
    groups: Mapping[str, Sequence[str]],
    tol: float,
    inclusion_cutoff: Optional[float] = None,
) -> MilpSolution:
    """Read a mixed-integer solution off both its indicators and its flux.

    Parameters
    ----------
    model : cobra.Model
        The model that was just optimized.
    directions : dict of {str: int}
        The oriented direction of each reaction.
    indicators : dict of {str: optlang.Variable}
        The binary inclusion variable of each indicator group.
    groups : dict of {str: sequence of str}
        The reactions each binary governs.
    tol : float
        The flux threshold, used to derive the inclusion cutoff.
    inclusion_cutoff : float, optional
        The absolute flux at which a reaction counts as used.

    Returns
    -------
    MilpSolution
        The reactions to keep, and the ones the indicators selected.

    """
    included = _included_from_indicators(directions, indicators, groups)
    carrying = _included_reactions(model, tol, inclusion_cutoff)
    selected = {key for key, var in indicators.items() if var.primal > 0.5}
    return MilpSolution(included=included | carrying, selected=selected)


def _free_reaction_ids(model: "Model", directions: Dict[str, int]) -> List[str]:
    """Return the identifiers of the reactions free to be dropped.

    Parameters
    ----------
    model : cobra.Model
        The model to inspect.
    directions : dict of {str: int}
        The oriented direction of each reaction.

    Returns
    -------
    list of str
        The identifiers of the reactions whose direction is 0, in model order.

    """
    return [rxn.id for rxn in model.reactions if directions.get(rxn.id, 0) == 0]


def _resolve_indicator_groups(
    free_ids: Sequence[str],
    indicator_reactions: Optional[Iterable[str]],
    indicator_groups: Optional[Mapping[str, Iterable[str]]],
) -> Tuple[Dict[str, List[str]], Set[str]]:
    """Work out which reactions each binary governs.

    Parameters
    ----------
    free_ids : sequence of str
        The reactions free to be dropped, in model order.
    indicator_reactions : iterable of str, optional
        Restrict the binaries to these reactions, one each.
    indicator_groups : dict of {str: iterable of str}, optional
        Give the reactions of each group a single shared binary.

    Returns
    -------
    dict of {str: list of str}
        The reactions each binary governs, keyed by group.
    set of str
        The free reactions with no binary, which are kept regardless.

    Raises
    ------
    SpectraError
        If both arguments are given, a group is empty, a reaction appears
        in two groups, or a named reaction is not free.

    Notes
    -----
    With neither argument every free reaction gets its own binary and the
    key is the reaction's own identifier, which is why `selected` reads as
    a set of reactions in that case.

    """
    if indicator_reactions is not None and indicator_groups is not None:
        raise SpectraError(
            "Pass indicator_reactions or indicator_groups, not both: one "
            "gives each named reaction its own binary and the other shares "
            "one binary between several, so giving both leaves it ambiguous "
            "which applies."
        )

    order = {rxn_id: index for index, rxn_id in enumerate(free_ids)}
    if indicator_groups is not None:
        groups: Dict[str, List[str]] = {}
        seen: Dict[str, str] = {}
        for key, members in indicator_groups.items():
            wanted = list(dict.fromkeys(members))
            if not wanted:
                raise SpectraError(
                    f"Indicator group {key!r} is empty, so its binary would "
                    f"govern nothing while still being counted."
                )
            unknown = [r for r in wanted if r not in order]
            if unknown:
                raise SpectraError(
                    f"indicator_groups must name free reactions of the model, "
                    f"but group {key!r} names {len(unknown)} that are not, for "
                    f"example {sorted(unknown)[:5]}. A reaction given a "
                    f"direction is already forced in and cannot also be "
                    f"selected over."
                )
            for rxn_id in wanted:
                if rxn_id in seen:
                    raise SpectraError(
                        f"{rxn_id!r} is in both group {seen[rxn_id]!r} and "
                        f"group {key!r}, so which binary governs it is "
                        f"ambiguous. Every reaction belongs to at most one."
                    )
                seen[rxn_id] = key
            groups[key] = sorted(wanted, key=order.__getitem__)
        return groups, set(free_ids) - set(seen)

    if indicator_reactions is None:
        return {rxn_id: [rxn_id] for rxn_id in free_ids}, set()

    wanted = set(indicator_reactions)
    unknown = wanted - set(free_ids)
    if unknown:
        raise SpectraError(
            f"indicator_reactions must be free reactions of the model, "
            f"but {len(unknown)} are not, for example "
            f"{sorted(unknown)[:5]}. A reaction given a direction is "
            f"already forced in and cannot also be selected over."
        )
    chosen = [rxn_id for rxn_id in free_ids if rxn_id in wanted]
    return {rxn_id: [rxn_id] for rxn_id in chosen}, set(free_ids) - wanted


def _add_absolute_value_vars(
    model: "Model", rxn_ids: List[str], upper_bound: Optional[float] = None
) -> Dict[str, "object"]:
    """Add a variable bounding the absolute flux of each given reaction.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on.
    rxn_ids : list of str
        The reactions to add an absolute-value variable for.
    upper_bound : float, optional
        An upper bound for the variables (default None, i.e. unbounded).

    Returns
    -------
    dict of {str: optlang.Variable}
        The variable bounding each reaction's absolute flux, such that
        ``t_i >= |v_i|``.

    """
    prob = model.problem
    variables = {}
    to_add = []

    for rxn_id in rxn_ids:
        reaction = model.reactions.get_by_id(rxn_id)
        var = prob.Variable(f"spectra_abs_{rxn_id}", lb=0.0, ub=upper_bound)
        variables[rxn_id] = var
        to_add.extend(
            [
                var,
                # t_i >= v_i and t_i >= -v_i, i.e. t_i >= |v_i|.
                prob.Constraint(
                    reaction.flux_expression - var,
                    name=f"spectra_abs_upper_{rxn_id}",
                    ub=0.0,
                ),
                prob.Constraint(
                    reaction.flux_expression + var,
                    name=f"spectra_abs_lower_{rxn_id}",
                    lb=0.0,
                ),
            ]
        )

    model.add_cons_vars(to_add)
    return variables


def _add_exclusion_constraints(
    model: "Model",
    indicators: Dict[str, "object"],
    previous_solutions: Optional[List[Set[str]]],
) -> None:
    """Forbid each previously found reaction set from recurring exactly.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on.
    indicators : dict of {str: optlang.Variable}
        The binary inclusion variable of each free reaction.
    previous_solutions : list of set of str, optional
        The reaction sets to exclude. Only their free reactions matter, since
        the directed ones are in every solution by construction.

    Notes
    -----
    For each previous solution, ``sum(z_i) <= |S| - 1`` over the indicators
    it switched on. That rules out the set itself **and every superset of
    it**, since a superset has all of ``S`` on and so hits the same bound;
    proper subsets stay available.

    Excluding supersets costs nothing here and is arguably what you want.
    Solutions of equal size can never contain one another, so no optimum is
    lost; what is pruned is only a larger answer that keeps everything an
    earlier one did.

    The sum has to be written over the indicators rather than over the
    reactions finally kept: a reaction swept in by the inclusion cutoff
    while its binary was off would otherwise appear in ``S`` without being
    counted, the sum would stay under the bound, and the same solution could
    recur. With `indicator_groups` the keys are groups rather than
    reactions, which is why :class:`MilpSolution` reports `selected` in
    those terms.

    """
    if not previous_solutions:
        return

    prob = model.problem
    constraints = []
    for index, solution in enumerate(previous_solutions):
        in_solution = [indicators[r] for r in solution if r in indicators]
        if not in_solution:
            continue
        constraints.append(
            prob.Constraint(
                add(in_solution),
                name=f"spectra_exclude_{index}",
                ub=len(in_solution) - 1,
            )
        )

    model.add_cons_vars(constraints)


def min_net_lp(
    model: "Model",
    directions: Dict[str, int],
    weights: Dict[str, float],
    tol: float,
    steady_state: bool = True,
    inclusion_cutoff: Optional[float] = None,
) -> Set[str]:
    """Extract a model by minimising the weighted total absolute flux.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from.
    directions : dict of {str: int}
        The oriented direction each reaction must carry flux in.
    weights : dict of {str: float}
        The weight of each reaction. Higher weights make a reaction less
        likely to be kept; they should be non-negative.
    tol : float
        The minimum absolute flux a directed reaction has to carry.
    steady_state : bool, optional
        Whether to enforce ``S v = 0`` rather than ``S v >= 0`` (default True).
    inclusion_cutoff : float, optional
        The absolute flux at which a reaction counts as part of the extracted
        model (default ``tol * 1e-7``). Raising it gives a smaller model but
        risks dropping reactions that are load-bearing for a mass balance; see
        :func:`_included_reactions`.

    Returns
    -------
    set of str
        The identifiers of the reactions in the extracted model.

    Notes
    -----
    A single LP minimising
    :math:`\\sum_i w_i t_i` subject to :math:`t_i \\ge |v_i|`, so it is an
    L1 relaxation of minimising the reaction count. Being a pure LP it is the
    cheapest formulation, at the cost of favouring many small fluxes over few
    large ones.

    """
    free_ids = _free_reaction_ids(model, directions)
    signs = reaction_signs(model)

    with model, relaxed_mass_balance(model, steady_state):
        apply_direction_bounds(model, directions, signs, tol)
        abs_vars = _add_absolute_value_vars(model, free_ids)

        model.objective = model.problem.Objective(Zero, direction="min")
        model.objective.set_linear_coefficients(
            {abs_vars[r]: float(weights[r]) for r in free_ids}
        )

        model.slim_optimize()
        _check_status(model, "minNetLP problem")
        return _included_reactions(model, tol, inclusion_cutoff)


def min_net_milp(
    model: "Model",
    directions: Dict[str, int],
    weights: Dict[str, float],
    tol: float,
    steady_state: bool = True,
    time_limit: Optional[float] = 7200.0,
    previous_solutions: Optional[List[Set[str]]] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
    inclusion_cutoff: Optional[float] = None,
    indicator_groups: Optional[Mapping[str, Iterable[str]]] = None,
) -> MilpSolution:
    """Extract a model by minimising the weighted reaction count.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from.
    directions : dict of {str: int}
        The oriented direction each reaction must carry flux in.
    weights : dict of {str: float}
        The weight of each reaction. Higher weights make a reaction less
        likely to be kept; they should be non-negative.
    tol : float
        The minimum absolute flux a directed reaction has to carry.
    steady_state : bool, optional
        Whether to enforce ``S v = 0`` rather than ``S v >= 0`` (default True).
    time_limit : float, optional
        The maximum time to spend solving, in seconds (default 7200). A
        solution found before the limit is used even if not proven optimal.
    previous_solutions : list of set of str, optional
        Reaction sets to exclude, so that a different model is returned
        (default None).
    indicator_reactions : iterable of str, optional
        Restrict the binaries to these reactions; the rest stay continuous
        and are kept whatever they do (default None, meaning every free
        reaction gets one). Choosing a subset is what makes a selection
        problem over a handful of reactions tractable on a large model: one
        binary per organism in a community, say, rather than one per
        reaction.
    inclusion_cutoff : float, optional
        The absolute flux at which a reaction without an indicator counts as
        used (default ``tol * 1e-7``).
    indicator_groups : dict of {str: iterable of str}, optional
        Give the reactions of each group **one shared binary** instead of
        one each, so that they are kept or dropped together and cost the
        group's weight once however many of them there are. Mutually
        exclusive with `indicator_reactions`. Any free reaction in no group
        is kept regardless, as with `indicator_reactions`.

    Returns
    -------
    MilpSolution
        The reactions to keep, and the indicator groups the solver
        switched on.

    Notes
    -----
    Each indicator gets a binary :math:`z_g` tied to the flux of every
    reaction it governs by :math:`lb_i z_g \\le v_i \\le ub_i z_g`, so
    :math:`z_g = 0` pins all of them to zero at once, and the objective
    minimises :math:`\\sum_g w_g z_g`. With the default one binary per
    reaction each group has a single member and this is the usual
    formulation. Unlike :func:`min_net_lp` it minimises the count exactly
    rather than its L1 relaxation, but it is a MILP and so far more
    expensive.

    Sharing a binary is what expresses a requirement that holds in several
    conditions at once. Replicate the network once per condition, group
    each reaction's copies together, and the objective counts reactions
    while the constraints are satisfied separately in every condition.
    That is what :func:`~spectra_cobra.minimal_reactome` does.

    Reactions left without an indicator keep their own bounds and are always
    included: nothing is deciding whether to keep them, so nothing should
    discard them.

    """
    free_ids = _free_reaction_ids(model, directions)
    groups, always_keep = _resolve_indicator_groups(
        free_ids, indicator_reactions, indicator_groups
    )
    missing = sorted(key for key in groups if key not in weights)
    if missing:
        raise SpectraError(
            f"{len(missing)} indicator group(s) have no weight, for example "
            f"{missing[:5]}. What a group costs to keep is looked up by its "
            f"own key, which for the default one-binary-per-reaction case is "
            f"the reaction identifier."
        )
    signs = reaction_signs(model)
    prob = model.problem

    with model, relaxed_mass_balance(model, steady_state):
        apply_direction_bounds(model, directions, signs, tol)

        indicators = {}
        to_add = []
        for key, members in groups.items():
            var = prob.Variable(f"spectra_z_{key}", type="binary")
            indicators[key] = var
            to_add.append(var)
            for rxn_id in members:
                reaction = model.reactions.get_by_id(rxn_id)
                # lb_i z_g <= v_i <= ub_i z_g. This pair is unaffected by the
                # orientation: negating both the flux and the swapped bounds
                # reproduces the same inequalities. One binary over several
                # reactions pins all of them at once.
                to_add.extend(
                    [
                        prob.Constraint(
                            reaction.flux_expression - reaction.lower_bound * var,
                            name=f"spectra_z_lower_{rxn_id}",
                            lb=0.0,
                        ),
                        prob.Constraint(
                            reaction.flux_expression - reaction.upper_bound * var,
                            name=f"spectra_z_upper_{rxn_id}",
                            ub=0.0,
                        ),
                    ]
                )
        model.add_cons_vars(to_add)
        _add_exclusion_constraints(model, indicators, previous_solutions)

        model.objective = prob.Objective(Zero, direction="min")
        model.objective.set_linear_coefficients(
            {indicators[key]: float(weights[key]) for key in groups}
        )
        _set_time_limit(model, time_limit)

        model.slim_optimize()
        _check_status(model, "minNetMILP problem", milp=True)
        solution = _milp_solution(
            model, directions, indicators, groups, tol, inclusion_cutoff
        )
        return MilpSolution(
            included=solution.included | always_keep, selected=solution.selected
        )


def trade_off(
    model: "Model",
    directions: Dict[str, int],
    weights: Dict[str, float],
    tol: float,
    steady_state: bool = True,
    time_limit: Optional[float] = 7200.0,
    previous_solutions: Optional[List[Set[str]]] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
    inclusion_cutoff: Optional[float] = None,
) -> MilpSolution:
    """Extract a model by maximising the weighted number of reactions kept.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from.
    directions : dict of {str: int}
        The oriented direction each reaction must carry flux in.
    weights : dict of {str: float}
        The weight of each reaction. Unlike the ``min_net`` formulations these
        may be **any real number**: a positive weight pushes a reaction into
        the model and a negative one pushes it out, so the solver trades
        positive evidence against negative.
    tol : float
        The minimum absolute flux a kept reaction has to carry.
    steady_state : bool, optional
        Whether to enforce ``S v = 0`` rather than ``S v >= 0`` (default True).
    time_limit : float, optional
        The maximum time to spend solving, in seconds (default 7200).
    previous_solutions : list of set of str, optional
        Reaction sets to exclude, so that a different model is returned
        (default None).
    indicator_reactions : iterable of str, optional
        Restrict the binaries to these reactions; the rest stay continuous
        and are kept whatever they do (default None, meaning every free
        reaction gets one).
    inclusion_cutoff : float, optional
        The absolute flux at which a reaction without an indicator counts as
        used (default ``tol * 1e-7``).

    Returns
    -------
    MilpSolution
        The reactions to keep, and the ones the indicators selected.

    Notes
    -----
    A free irreversible reaction gets one binary
    :math:`z_i` with :math:`\\varepsilon z_i \\le \\hat v_i \\le ub_i z_i`,
    where :math:`\\hat v` is the oriented flux, so being included means
    carrying at least :math:`\\varepsilon`. A free *reversible* reaction gets
    two more binaries, :math:`a_i` for running forward and :math:`b_i` for
    running backward, constrained by :math:`a_i + b_i = z_i` so that at most
    one direction is active:

    .. math::

        \\hat v_i &\\ge \\varepsilon a_i + lb_i b_i \\\\
        \\hat v_i &\\le ub_i a_i - \\varepsilon b_i

    With :math:`a_i = b_i = 0` these pin the flux to zero.

    """
    free_ids = _free_reaction_ids(model, directions)
    # tradeOff gives a reversible reaction two further binaries of its own,
    # so a shared one would have to govern those too; it takes
    # indicator_reactions only, and the resolver is used just for the
    # validation and the always-keep set.
    one_each, always_keep = _resolve_indicator_groups(
        free_ids, indicator_reactions, None
    )
    indicator_ids = list(one_each)
    signs = reaction_signs(model)
    prob = model.problem

    with model, relaxed_mass_balance(model, steady_state):
        apply_direction_bounds(model, directions, signs, tol)

        indicators = {}
        to_add = []
        for rxn_id in indicator_ids:
            reaction = model.reactions.get_by_id(rxn_id)
            sign = signs[rxn_id]
            lower, upper = oriented_bounds(reaction, sign)
            oriented_flux = sign * reaction.flux_expression

            included = prob.Variable(f"spectra_z_{rxn_id}", type="binary")
            indicators[rxn_id] = included
            to_add.append(included)

            if not is_reversible(reaction, sign):
                # tol * z <= oriented flux <= ub * z
                to_add.extend(
                    [
                        prob.Constraint(
                            oriented_flux - tol * included,
                            name=f"spectra_irr_lower_{rxn_id}",
                            lb=0.0,
                        ),
                        prob.Constraint(
                            oriented_flux - upper * included,
                            name=f"spectra_irr_upper_{rxn_id}",
                            ub=0.0,
                        ),
                    ]
                )
                continue

            forward_on = prob.Variable(f"spectra_a_{rxn_id}", type="binary")
            reverse_on = prob.Variable(f"spectra_b_{rxn_id}", type="binary")
            to_add.extend(
                [
                    forward_on,
                    reverse_on,
                    # a + b = z: at most one direction may be active.
                    prob.Constraint(
                        forward_on + reverse_on - included,
                        name=f"spectra_split_{rxn_id}",
                        lb=0.0,
                        ub=0.0,
                    ),
                    prob.Constraint(
                        oriented_flux - tol * forward_on - lower * reverse_on,
                        name=f"spectra_rev_lower_{rxn_id}",
                        lb=0.0,
                    ),
                    prob.Constraint(
                        oriented_flux - upper * forward_on + tol * reverse_on,
                        name=f"spectra_rev_upper_{rxn_id}",
                        ub=0.0,
                    ),
                ]
            )

        model.add_cons_vars(to_add)
        _add_exclusion_constraints(model, indicators, previous_solutions)

        model.objective = prob.Objective(Zero, direction="max")
        model.objective.set_linear_coefficients(
            {indicators[r]: float(weights[r]) for r in indicator_ids}
        )
        _set_time_limit(model, time_limit)

        model.slim_optimize()
        _check_status(model, "tradeOff problem", milp=True)
        solution = _milp_solution(
            model, directions, indicators, one_each, tol, inclusion_cutoff
        )
        return MilpSolution(
            included=solution.included | always_keep, selected=solution.selected
        )


def growth_optim(
    model: "Model",
    directions: Dict[str, int],
    weights: Dict[str, float],
    tol: float,
    steady_state: bool = True,
    inclusion_cutoff: Optional[float] = None,
) -> Set[str]:
    """Extract a model by maximising growth while penalising total flux.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from. Its objective identifies the biomass
        reactions.
    directions : dict of {str: int}
        The oriented direction each reaction must carry flux in.
    weights : dict of {str: float}
        The weight of each reaction. For the objective reactions the weight is
        the reward on their flux; for every other reaction it is the penalty
        on its absolute flux. They should be non-negative.
    tol : float
        The minimum absolute flux a directed reaction has to carry.
    steady_state : bool, optional
        Whether to enforce ``S v = 0`` rather than ``S v >= 0`` (default True).
    inclusion_cutoff : float, optional
        The absolute flux at which a reaction counts as part of the extracted
        model (default ``tol * 1e-7``). Raising it gives a smaller model but
        risks dropping reactions that are load-bearing for a mass balance; see
        :func:`_included_reactions`.

    Returns
    -------
    set of str
        The identifiers of the reactions in the extracted model.

    Raises
    ------
    SpectraError
        If the model has no objective reaction to grow.

    Notes
    -----
    Maximises
    :math:`\\sum_{i \\in c} w_i v_i - \\sum_{i \\notin c} w_i t_i` with
    :math:`t_i \\ge |v_i|`, where :math:`c` are the objective reactions. Note
    that the reward on the objective reactions is their *weight*, not their
    coefficient in ``model.objective``: the objective is used only to find
    which reactions those are.

    Unlike the other formulations this one penalises the absolute flux of
    every non-objective reaction, core reactions included, not just the free
    ones.

    """
    objective_ids = {
        rxn.id for rxn, coef in linear_reaction_coefficients(model).items() if coef
    }
    if not objective_ids:
        raise SpectraError(
            "growth_optim needs an objective reaction to maximise, but this "
            "model's objective has no nonzero linear coefficients."
        )

    penalised_ids = [rxn.id for rxn in model.reactions if rxn.id not in objective_ids]
    signs = reaction_signs(model)

    with model, relaxed_mass_balance(model, steady_state):
        apply_direction_bounds(model, directions, signs, tol)
        abs_vars = _add_absolute_value_vars(model, penalised_ids)

        model.objective = model.problem.Objective(Zero, direction="max")
        coefficients = {abs_vars[r]: -float(weights[r]) for r in penalised_ids}
        for rxn_id in objective_ids:
            reaction = model.reactions.get_by_id(rxn_id)
            coefficients[reaction.forward_variable] = float(weights[rxn_id])
            coefficients[reaction.reverse_variable] = -float(weights[rxn_id])
        model.objective.set_linear_coefficients(coefficients)

        model.slim_optimize()
        _check_status(model, "growthOptim problem")
        return _included_reactions(model, tol, inclusion_cutoff)


def _set_time_limit(model: "Model", time_limit: Optional[float]) -> None:
    """Ask the solver to stop after the given time, if it supports that.

    Parameters
    ----------
    model : cobra.Model
        The model whose solver to configure.
    time_limit : float, optional
        The limit in seconds, or None to leave the solver's default in place.

    Notes
    -----
    Not every optlang interface exposes a timeout, so a solver that does not
    is logged and left alone rather than treated as an error.

    """
    if time_limit is None:
        return
    try:
        # glpk's setter is typed as an int number of seconds, so round up
        # rather than handing it a float it will reject.
        model.solver.configuration.timeout = int(math.ceil(time_limit))
    except (AttributeError, TypeError, ValueError) as error:
        logger.warning(
            "Could not set a time limit of %.4gs on the %s solver: %s",
            time_limit,
            model.problem.__name__,
            error,
        )
