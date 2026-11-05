"""The three auxiliary LPs that SPECTRA's iterative phases are built from.

All three share a shape: an auxiliary variable per reaction of interest is
tied to that reaction's oriented flux by a one-sided constraint, and a
randomly weighted sum of those auxiliary variables is pushed as far as the
flux bounds allow. The randomised weights break ties between reactions, so
repeated calls explore different corners of the optimal face.

The three differ only in which reactions they cover and which way they push:

================  ==========================  ==========================
Routine           Covers                      Pushes oriented flux
================  ==========================  ==========================
``forward_cc``    every reaction given        up, towards ``+tol``
``forward``       reversible reactions given  up, towards ``+tol``
``reverse``       reversible reactions given  down, towards ``-tol``
================  ==========================  ==========================

``forward`` additionally pins the *irreversible* reactions it is given to at
least ``tol`` through their bounds rather than through an auxiliary variable,
which is what makes it able to fail: an infeasible core set has no such flux
distribution.
"""

from logging import getLogger
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Set

import numpy as np
from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

from ._orientation import is_reversible, oriented_bounds
from .exceptions import SpectraInfeasibleCoreError, SpectraSolverError

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)

#: Lower and upper limit of the uniformly sampled objective coefficients.
WEIGHT_RANGE = (1.0, 1.1)


def _fluxes(model: "Model") -> Dict[str, float]:
    """Return the current flux of every reaction, keyed by identifier.

    Parameters
    ----------
    model : cobra.Model
        A model whose solver holds an optimal solution.

    Returns
    -------
    dict of {str: float}
        The flux through each reaction.

    """
    return {rxn.id: rxn.flux for rxn in model.reactions}


def _zero_fluxes(model: "Model") -> Dict[str, float]:
    """Return a flux of zero for every reaction, keyed by identifier.

    Parameters
    ----------
    model : cobra.Model
        The model to describe.

    Returns
    -------
    dict of {str: float}
        A zero flux for each reaction.

    """
    return {rxn.id: 0.0 for rxn in model.reactions}


def _solve(
    model: "Model",
    description: str,
    raise_on_failure: bool,
) -> Optional[Dict[str, float]]:
    """Optimize the model and return its fluxes, or None if it failed.

    Parameters
    ----------
    model : cobra.Model
        The model to optimize.
    description : str
        What is being solved, used in the warning or error message.
    raise_on_failure : bool
        Whether a non-optimal status is an error rather than a warning.

    Returns
    -------
    dict of {str: float} or None
        The fluxes if the solver reached optimality, None otherwise.

    Raises
    ------
    SpectraSolverError
        If the solver did not reach optimality and `raise_on_failure` is set.

    """
    model.slim_optimize()
    status = model.solver.status
    if status == OPTIMAL:
        return _fluxes(model)

    if raise_on_failure:
        raise SpectraSolverError(
            f"The {description} terminated with status {status!r} instead of "
            f"reaching optimality.",
            status=status,
        )
    logger.warning(
        "The %s terminated with status %r; treating its reactions as carrying "
        "no flux.",
        description,
        status,
    )
    return None


def _push(
    model: "Model",
    rxn_ids: Iterable[str],
    signs: Dict[str, float],
    tol: float,
    rng: np.random.Generator,
    reverse: bool,
    description: str,
    raise_on_failure: bool = False,
) -> Optional[Dict[str, float]]:
    """Push the oriented flux of the given reactions towards the threshold.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on. Call this inside a ``with model:`` block.
    rxn_ids : iterable of str
        The reactions to tie an auxiliary variable to.
    signs : dict of {str: float}
        The orientation of each reaction.
    tol : float
        The flux threshold the auxiliary variables are bounded by.
    rng : numpy.random.Generator
        The source of the objective coefficients.
    reverse : bool
        Whether to push the oriented flux down rather than up.
    description : str
        What is being solved, for the warning or error message.
    raise_on_failure : bool, optional
        Whether a non-optimal status is an error (default False).

    Returns
    -------
    dict of {str: float} or None
        The resulting fluxes, or None if the solver did not reach optimality
        and `raise_on_failure` is not set.

    """
    # Sorted, not merely listed: the objective coefficients below are drawn
    # in one block and zipped against the variables positionally, so the
    # order decides which reaction gets which weight. Callers pass a set of
    # identifiers, whose iteration order varies with PYTHONHASHSEED, and that
    # would make a seeded run irreproducible between processes.
    rxn_ids = sorted(rxn_ids)
    prob = model.problem
    aux_vars = []
    constraints = []

    # The constraint is sign_i * v_i - z_i, bounded on one side. Writing that
    # as an expression makes optlang build it symbolically, and sympy then
    # dominates the run: 10600 of them cost 19.1 s to construct against 2.3 s
    # for the same constraints declared empty and filled in by coefficient.
    # ``v_i`` is itself ``forward - reverse``, so the coefficients are
    # sign_i on the forward variable, -sign_i on the reverse one, and -1 on
    # the auxiliary variable.
    coefficients = []
    for rxn_id in rxn_ids:
        reaction = model.reactions.get_by_id(rxn_id)
        sign = signs[rxn_id]

        if reverse:
            # z_i >= -tol and v_i <= z_i, minimised: drives v_i down to -tol.
            # Unbounded above on purpose: a finite bound would silently
            # bind on a model with large flux bounds.
            var = prob.Variable(f"spectra_aux_{rxn_id}", lb=-tol, ub=None)
            constraint = prob.Constraint(
                Zero, name=f"spectra_aux_cons_{rxn_id}", ub=0.0
            )
        else:
            # z_i <= tol and v_i >= z_i, maximised: drives v_i up to tol.
            # Unbounded below on purpose, so that a reaction unable to go
            # positive merely fails to contribute rather than making the
            # whole LP infeasible.
            var = prob.Variable(f"spectra_aux_{rxn_id}", lb=None, ub=tol)
            constraint = prob.Constraint(
                Zero, name=f"spectra_aux_cons_{rxn_id}", lb=0.0
            )
        constraints.append(constraint)
        aux_vars.append(var)
        coefficients.append(
            (
                constraint,
                {
                    reaction.forward_variable: sign,
                    reaction.reverse_variable: -sign,
                    var: -1.0,
                },
            )
        )

    to_add = aux_vars + constraints
    model.add_cons_vars(to_add)
    # The coefficients can only be set once the constraints belong to a
    # problem.
    model.solver.update()
    for constraint, terms in coefficients:
        constraint.set_linear_coefficients(terms)
    # Reversed so the context rollback deletes the tail of optlang's variable
    # container first: removing from the front re-indexes everything after it,
    # which is quadratic in the number of auxiliary variables.
    to_add.reverse()
    model.objective = prob.Objective(Zero, direction="min" if reverse else "max")
    if aux_vars:
        weights = rng.uniform(*WEIGHT_RANGE, size=len(aux_vars))
        model.objective.set_linear_coefficients(dict(zip(aux_vars, weights)))

    # The LP is solved even with no auxiliary variables, because `forward`
    # pins its irreversible reactions through their bounds rather than through
    # the objective: the flux it is after is there regardless.
    return _solve(model, description, raise_on_failure)


def forward_cc(
    model: "Model",
    rxn_ids: Iterable[str],
    signs: Dict[str, float],
    tol: float,
    rng: np.random.Generator,
) -> Optional[Dict[str, float]]:
    """Drive as many of the given reactions as possible to positive flux.

    Every reaction given gets an auxiliary variable, and none of them is
    forced to carry flux, so the LP is always feasible as long as the model
    itself is.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on. Call this inside a ``with model:`` block.
    rxn_ids : iterable of str
        The reactions to drive.
    signs : dict of {str: float}
        The orientation of each reaction.
    tol : float
        The flux threshold.
    rng : numpy.random.Generator
        The source of the objective coefficients.

    Returns
    -------
    dict of {str: float} or None
        The resulting fluxes, or None if the solver did not reach optimality.

    """
    return _push(
        model,
        rxn_ids,
        signs,
        tol,
        rng,
        reverse=False,
        description="forward consistency LP",
    )


def reverse(
    model: "Model",
    rxn_ids: Iterable[str],
    signs: Dict[str, float],
    tol: float,
    rng: np.random.Generator,
    raise_on_failure: bool = False,
) -> Optional[Dict[str, float]]:
    """Drive as many of the given reversible reactions as possible negative.

    Only the reversible reactions among `rxn_ids` are considered, since an
    irreversible one cannot carry negative oriented flux by definition.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on. Call this inside a ``with model:`` block.
    rxn_ids : iterable of str
        The reactions to consider; the irreversible ones are skipped.
    signs : dict of {str: float}
        The orientation of each reaction.
    tol : float
        The flux threshold.
    rng : numpy.random.Generator
        The source of the objective coefficients.
    raise_on_failure : bool, optional
        Whether a non-optimal status is an error (default False).

    Returns
    -------
    dict of {str: float} or None
        The resulting fluxes, or None if the solver did not reach optimality
        and `raise_on_failure` is not set.

    """
    reversible = [
        rxn_id
        for rxn_id in rxn_ids
        if is_reversible(model.reactions.get_by_id(rxn_id), signs[rxn_id])
    ]
    return _push(
        model,
        reversible,
        signs,
        tol,
        rng,
        reverse=True,
        description="reverse LP",
        raise_on_failure=raise_on_failure,
    )


def forward(
    model: "Model",
    rxn_ids: Iterable[str],
    signs: Dict[str, float],
    tol: float,
    rng: np.random.Generator,
) -> Dict[str, float]:
    """Drive the given reactions to positive flux, pinning irreversible ones.

    This is the variant used during model extraction. It differs from
    :func:`forward_cc` in that the *irreversible* reactions given are required
    to carry at least `tol` through their lower bound, rather than merely
    rewarded for doing so. Only the reversible ones get an auxiliary
    variable.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on. Call this inside a ``with model:`` block.
    rxn_ids : iterable of str
        The reactions to drive. The irreversible ones among them are pinned.
    signs : dict of {str: float}
        The orientation of each reaction.
    tol : float
        The flux threshold.
    rng : numpy.random.Generator
        The source of the objective coefficients.

    Returns
    -------
    dict of {str: float}
        The resulting fluxes.

    Raises
    ------
    SpectraInfeasibleCoreError
        If no flux distribution lets every pinned reaction carry `tol` at
        once. Raising is deliberate: returning a ``NaN`` solution here leaves
        the caller looping forever.

    """
    rxn_ids = list(rxn_ids)
    reversible, irreversible = [], []
    for rxn_id in rxn_ids:
        reaction = model.reactions.get_by_id(rxn_id)
        if is_reversible(reaction, signs[rxn_id]):
            reversible.append(rxn_id)
        else:
            irreversible.append(rxn_id)

    # Pin the irreversible reactions through their bounds:
    #     lb(irr_core) = max(tol, lb(irr_core))
    for rxn_id in irreversible:
        reaction = model.reactions.get_by_id(rxn_id)
        sign = signs[rxn_id]
        lower, upper = oriented_bounds(reaction, sign)
        lower = max(tol, lower)
        if sign > 0:
            reaction.bounds = (lower, upper)
        else:
            reaction.bounds = (-upper, -lower)

    try:
        fluxes = _push(
            model,
            reversible,
            signs,
            tol,
            rng,
            reverse=False,
            description="forward extraction LP",
            raise_on_failure=True,
        )
    except SpectraSolverError as error:
        raise SpectraInfeasibleCoreError(
            f"No flux distribution lets all {len(irreversible)} pinned "
            f"irreversible core reactions carry at least tol={tol:.3g} at "
            f"once (solver status {error.status!r}). Either the core set is "
            f"mutually inconsistent, or some of its reactions are blocked in "
            f"this model; spectra_ccme handles the latter by dropping them."
        ) from error

    return fluxes


def carrying_flux(fluxes: Dict[str, float], cutoff: float) -> Set[str]:
    """Return the reactions whose absolute flux reaches the cutoff.

    Parameters
    ----------
    fluxes : dict of {str: float}
        The flux through each reaction.
    cutoff : float
        The threshold an absolute flux has to reach.

    Returns
    -------
    set of str
        The identifiers of the reactions at or above the cutoff.

    """
    return {rxn_id for rxn_id, flux in fluxes.items() if abs(flux) >= cutoff}
