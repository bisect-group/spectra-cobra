"""Sparsify a flux distribution, ``FluxReducer``."""

from logging import getLogger
from typing import TYPE_CHECKING, Dict

from optlang.symbolics import Zero

from ._orientation import (
    STOICHIOMETRY,
    is_reversible,
    reaction_signs,
    relaxed_mass_balance,
    validate_consistency_type,
)
from .formulations import _add_absolute_value_vars, _check_status

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)


def flux_reducer(
    model: "Model",
    consistency_type: str = STOICHIOMETRY,
) -> Dict[str, float]:
    """Find a sparse flux distribution for a model.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on.
    consistency_type : {"stoichiometry", "topology"}, optional
        Whether to assume a steady state (``S v = 0``) or an accumulation
        condition (``S v >= 0``) (default "stoichiometry").

    Returns
    -------
    dict of {str: float}
        The flux through each reaction.

    Notes
    -----
    This ports ``FluxReducer.m``, which minimises the sum of absolute fluxes.
    The objective is built so that an irreversible reaction is charged for its
    own flux directly, while a reversible one is charged through an
    absolute-value variable, which spares the problem an extra variable per
    irreversible reaction:

    .. math::

        \\text{minimize} \\quad \\sum_{i \\notin rev} \\hat v_i +
                                \\sum_{i \\in rev} t_i
        \\quad\\text{s.t.}\\quad t_i \\ge |v_i|

    The oriented flux :math:`\\hat v_i` is non-negative for an irreversible
    reaction, so the first sum is already an absolute value.

    """
    steady_state = validate_consistency_type(consistency_type)
    signs = reaction_signs(model)

    reversible_ids, irreversible_ids = [], []
    for rxn in model.reactions:
        if is_reversible(rxn, signs[rxn.id]):
            reversible_ids.append(rxn.id)
        else:
            irreversible_ids.append(rxn.id)

    with model, relaxed_mass_balance(model, steady_state):
        abs_vars = _add_absolute_value_vars(model, reversible_ids)

        model.objective = model.problem.Objective(Zero, direction="min")
        coefficients = {abs_vars[r]: 1.0 for r in reversible_ids}
        for rxn_id in irreversible_ids:
            reaction = model.reactions.get_by_id(rxn_id)
            sign = signs[rxn_id]
            coefficients[reaction.forward_variable] = sign
            coefficients[reaction.reverse_variable] = -sign
        model.objective.set_linear_coefficients(coefficients)

        model.slim_optimize()
        _check_status(model, "FluxReducer problem")
        return {rxn.id: rxn.flux for rxn in model.reactions}
