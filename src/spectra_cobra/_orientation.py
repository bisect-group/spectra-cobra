"""Shared machinery for orienting reactions and relaxing the mass balance.

Every routine needs a canonical "try this way first" direction for each
reaction, which classically is arranged by negating the stoichiometric
columns of reactions that can only carry negative flux and swapping their
bounds. Mutating a :class:`cobra.Model` that way would be both expensive and
surprising, so this module expresses the same thing as a per-reaction sign
applied to the reaction's flux expression. A reaction's
*oriented* flux is ``sign * flux``, and the two formulations are equivalent
because negating a column and negating its variable cancel out:

.. math::

    \\sum_i (\\pm S_{:,i})(\\pm v_i) = \\sum_i S_{:,i} v_i

so the mass balance is untouched by the orientation, and only the terms that
single out a direction need the sign.
"""

from contextlib import contextmanager
from typing import TYPE_CHECKING, Dict, Iterator, Optional, Tuple

from .exceptions import SpectraError

if TYPE_CHECKING:
    from cobra.core import Model, Reaction


STOICHIOMETRY = "stoichiometry"
TOPOLOGY = "topology"
CONSISTENCY_TYPES = (STOICHIOMETRY, TOPOLOGY)


def reaction_sign(reaction: "Reaction") -> float:
    """Return the orientation a reaction is examined in.

    Parameters
    ----------
    reaction : cobra.Reaction
        The reaction to orient.

    Returns
    -------
    float
        -1.0 if the reaction can only carry negative flux, 1.0 otherwise,
        i.e. the sign is negative exactly when the upper bound is at most
        zero.

    """
    return -1.0 if reaction.upper_bound <= 0.0 else 1.0


def reaction_signs(model: "Model") -> Dict[str, float]:
    """Return the orientation of every reaction, keyed by reaction identifier.

    Parameters
    ----------
    model : cobra.Model
        The model to orient.

    Returns
    -------
    dict of {str: float}
        The orientation of each reaction, as returned by :func:`reaction_sign`.

    """
    return {rxn.id: reaction_sign(rxn) for rxn in model.reactions}


def oriented_bounds(reaction: "Reaction", sign: Optional[float] = None) -> Tuple:
    """Return a reaction's bounds in its oriented frame.

    Parameters
    ----------
    reaction : cobra.Reaction
        The reaction whose bounds to orient.
    sign : float, optional
        The reaction's orientation. Computed with :func:`reaction_sign` when
        not given (default None).

    Returns
    -------
    tuple of (float, float)
        The lower and upper bound on ``sign * flux``. For an unflipped
        reaction these are simply its own bounds; for a flipped one they are
        its bounds negated and swapped.

    """
    if sign is None:
        sign = reaction_sign(reaction)
    scaled = (sign * reaction.lower_bound, sign * reaction.upper_bound)
    return min(scaled), max(scaled)


def is_reversible(reaction: "Reaction", sign: Optional[float] = None) -> bool:
    """Return whether a reaction can carry flux in its reverse orientation.

    Parameters
    ----------
    reaction : cobra.Reaction
        The reaction to inspect.
    sign : float, optional
        The reaction's orientation (default None, see :func:`reaction_sign`).

    Returns
    -------
    bool
        Whether the oriented lower bound is negative. Reversibility is
        judged *after* orienting, so a reaction restricted to negative flux
        counts as irreversible.

    """
    return oriented_bounds(reaction, sign)[0] < 0.0


def validate_consistency_type(consistency_type: str) -> bool:
    """Return whether to enforce a steady state, validating the given name.

    Parameters
    ----------
    consistency_type : {"stoichiometry", "topology"}
        Whether to assume a steady state (``S v = 0``) or an accumulation
        condition (``S v >= 0``).

    Returns
    -------
    bool
        True for ``"stoichiometry"`` and False for ``"topology"``.

    Raises
    ------
    SpectraError
        If `consistency_type` is neither of the two accepted values. Failing
        here is deliberate: an unrecognised value that is allowed through
        fails later and less helpfully.

    """
    if consistency_type not in CONSISTENCY_TYPES:
        raise SpectraError(
            f"consistency_type must be one of {CONSISTENCY_TYPES}, "
            f"not {consistency_type!r}."
        )
    return consistency_type == STOICHIOMETRY


@contextmanager
def relaxed_mass_balance(model: "Model", steady_state: bool) -> Iterator[None]:
    """Relax the mass balance to an accumulation condition, then restore it.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on.
    steady_state : bool
        When True nothing is changed and ``S v = 0`` stays in force. When
        False every metabolite constraint is relaxed to ``S v >= 0``, which
        is the topology mode.

    Yields
    ------
    None

    Notes
    -----
    A :class:`cobra.Model` context manager does not track edits made directly
    to the bounds of existing constraints, so the original bounds are captured
    and restored here instead.

    """
    if steady_state:
        yield
        return

    original = {}
    try:
        for metabolite in model.metabolites:
            constraint = model.constraints.get(metabolite.id)
            original[metabolite.id] = (constraint.lb, constraint.ub)
            # S v >= 0: keep the lower bound, drop the upper one.
            constraint.ub = None
        model.solver.update()
        yield
    finally:
        for met_id, (lb, ub) in original.items():
            constraint = model.constraints.get(met_id)
            constraint.lb = lb
            constraint.ub = ub
        model.solver.update()


def apply_direction_bounds(
    model: "Model",
    directions: Dict[str, int],
    signs: Dict[str, float],
    tol: float,
) -> None:
    """Force the reactions with a known direction to carry flux that way.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on. Call this inside a ``with model:`` block, so
        that the bound changes are reverted afterwards.
    directions : dict of {str: int}
        The oriented direction of each reaction: 1 to force positive oriented
        flux, -1 to force negative oriented flux, and 0 to leave it free.
    signs : dict of {str: float}
        The orientation of each reaction.
    tol : float
        The minimum absolute flux a directed reaction has to carry.

    Notes
    -----
    Every formulation shares this bound rewrite, in the oriented frame::

        lb(direction == 1) = max(tol, lb(direction == 1));
        ub(direction == -1) = -tol;

    expressed here in the model's own flux frame instead. Note the asymmetry:
    for ``direction == -1`` the upper bound is *assigned* ``-tol`` rather than
    taking a minimum, so a reaction whose oriented upper bound already sits
    below ``-tol`` has it raised.

    """
    for rxn_id, direction in directions.items():
        if direction == 0:
            continue

        reaction = model.reactions.get_by_id(rxn_id)
        sign = signs[rxn_id]
        lower, upper = oriented_bounds(reaction, sign)

        if direction == 1:
            lower = max(tol, lower)
        else:
            upper = -tol

        # Map the oriented bounds back onto the reaction's own frame.
        if sign > 0:
            reaction.bounds = (lower, upper)
        else:
            reaction.bounds = (-upper, -lower)
