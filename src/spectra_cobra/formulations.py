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
from logging import getLogger
from typing import TYPE_CHECKING, Dict, List, Optional, Set

from cobra.util.solver import linear_reaction_coefficients
from optlang.interface import FEASIBLE, OPTIMAL
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

#: Statuses a MILP may stop at and still carry a usable solution. The MATLAB
#: implementation accepts ``stat == 1 || stat == 3``, i.e. optimal or a
#: feasible solution found before the time limit.
ACCEPTABLE_MILP_STATUSES = (OPTIMAL, FEASIBLE)

#: The fraction of `tol` below which a flux counts as zero when reading the
#: extracted reaction set off a solution, matching ``tol * 1e-7`` in MATLAB.
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
        logger.warning(
            "The %s stopped at status %r rather than optimality; the returned "
            "model is feasible but may not be minimal.",
            description,
            status,
        )
    return status


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
        (default ``tol * 1e-7``, as in MATLAB).

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
    return {rxn.id for rxn in model.reactions if abs(rxn.flux) >= inclusion_cutoff}


def _included_from_indicators(
    directions: Dict[str, int], indicators: Dict[str, "object"]
) -> Set[str]:
    """Return the reactions a mixed-integer solution keeps.

    Parameters
    ----------
    directions : dict of {str: int}
        The oriented direction of each reaction. The directed ones are kept
        unconditionally, since they are forced to carry flux.
    indicators : dict of {str: optlang.Variable}
        The binary inclusion variable of each free reaction.

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
    kept.update(rxn_id for rxn_id, var in indicators.items() if var.primal > 0.5)
    return kept


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
    This ports the ``prevSols`` block shared by ``minNet.m`` and
    ``tradeOff.m``: for each previous solution, ``sum(z_i) <= |S| - 1`` over
    the free reactions it contained, which rules out that exact set while
    allowing any subset or superset.

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
    This is the ``minNetLP`` branch of ``minNet.m``: a single LP minimising
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
) -> Set[str]:
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

    Returns
    -------
    set of str
        The identifiers of the reactions in the extracted model.

    Notes
    -----
    This is the ``minNetMILP`` branch of ``minNet.m``. Each free reaction gets
    a binary variable :math:`z_i` tied to its flux by
    :math:`lb_i z_i \\le v_i \\le ub_i z_i`, so :math:`z_i = 0` pins the flux
    to zero, and the objective minimises :math:`\\sum_i w_i z_i`. Unlike
    :func:`min_net_lp` this minimises the count exactly rather than its L1
    relaxation, but it is a MILP and so far more expensive.

    """
    free_ids = _free_reaction_ids(model, directions)
    signs = reaction_signs(model)
    prob = model.problem

    with model, relaxed_mass_balance(model, steady_state):
        apply_direction_bounds(model, directions, signs, tol)

        indicators = {}
        to_add = []
        for rxn_id in free_ids:
            reaction = model.reactions.get_by_id(rxn_id)
            var = prob.Variable(f"spectra_z_{rxn_id}", type="binary")
            indicators[rxn_id] = var
            # lb_i z_i <= v_i <= ub_i z_i. This pair is unaffected by the
            # orientation: negating both the flux and the swapped bounds
            # reproduces the same inequalities.
            to_add.extend(
                [
                    var,
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
            {indicators[r]: float(weights[r]) for r in free_ids}
        )
        _set_time_limit(model, time_limit)

        model.slim_optimize()
        _check_status(model, "minNetMILP problem", milp=True)
        return _included_from_indicators(directions, indicators)


def trade_off(
    model: "Model",
    directions: Dict[str, int],
    weights: Dict[str, float],
    tol: float,
    steady_state: bool = True,
    time_limit: Optional[float] = 7200.0,
    previous_solutions: Optional[List[Set[str]]] = None,
) -> Set[str]:
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

    Returns
    -------
    set of str
        The identifiers of the reactions in the extracted model.

    Notes
    -----
    This ports ``tradeOff.m``. A free irreversible reaction gets one binary
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
    signs = reaction_signs(model)
    prob = model.problem

    with model, relaxed_mass_balance(model, steady_state):
        apply_direction_bounds(model, directions, signs, tol)

        indicators = {}
        to_add = []
        for rxn_id in free_ids:
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
            {indicators[r]: float(weights[r]) for r in free_ids}
        )
        _set_time_limit(model, time_limit)

        model.slim_optimize()
        _check_status(model, "tradeOff problem", milp=True)
        return _included_from_indicators(directions, indicators)


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
    This ports ``growthOptim.m``, which maximises
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
