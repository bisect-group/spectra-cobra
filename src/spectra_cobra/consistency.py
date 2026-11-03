"""The SPECTRA consistency check."""

from logging import getLogger
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Tuple

import numpy as np
from optlang.symbolics import Zero

from ._copy import subset_model
from ._lp import (
    _solve,
    build_push_block,
    carrying_flux,
    draw_weights,
    retarget_objective,
)
from ._orientation import (
    STOICHIOMETRY,
    is_reversible,
    reaction_signs,
    relaxed_mass_balance,
    validate_consistency_type,
)
from .formulations import DEFAULT_SOLVER_TOLERANCE, _apply_tolerances

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)


def consistent_reaction_ids(
    model: "Model",
    tol: float = 1e-4,
    consistency_type: str = STOICHIOMETRY,
    detection_cutoff: Optional[float] = None,
    seed: Optional[int] = None,
    feasibility_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    optimality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
) -> Tuple[Set[str], int]:
    """Find the flux consistent reactions of a model.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on.
    tol : float, optional
        The flux threshold, i.e. the magnitude the LPs drive the reaction
        fluxes to (default 1e-4).
    consistency_type : {"stoichiometry", "topology"}, optional
        Whether to assume a steady state (``S v = 0``) or an accumulation
        condition (``S v >= 0``) (default "stoichiometry").
    detection_cutoff : float, optional
        The absolute flux at which a reaction counts as carrying flux. The
        default requires a reaction to reach ``tol`` itself; pass
        ``model.tolerance`` instead to count any nonzero flux, matching
        :func:`cobra.flux_analysis.fastcc` (default ``0.99 * tol``).
    seed : int, optional
        A seed for the random objective coefficients, making the result
        reproducible (default None).
    feasibility_tolerance : float, optional
        How far the mass balance may be violated and still count as
        satisfied. A reaction counts as consistent when an LP drives it to
        ``tol``, so a feasibility tolerance near ``tol`` lets a reaction pass
        on residual rather than on real flux (default 1e-9, the lowest most
        solvers accept). Pass None to leave the solver's own in place.
    optimality_tolerance : float, optional
        How far a reduced cost may be from optimal and still count as optimal
        (default 1e-9). None leaves the solver's own in place.

    Returns
    -------
    tuple of (set of str, int)
        The identifiers of the consistent reactions and the number of LPs
        that were solved to find them.

    """
    steady_state = validate_consistency_type(consistency_type)
    if detection_cutoff is None:
        detection_cutoff = 0.99 * tol
    rng = np.random.default_rng(seed)

    # The consistency check has no binaries, so there is no integrality
    # tolerance to set.
    _apply_tolerances(model, feasibility_tolerance, optimality_tolerance)
    signs = reaction_signs(model)
    # Reactions that have not been shown to carry flux yet.
    rxns_to_check = set(signs)
    n_lps = 0

    with model, relaxed_mass_balance(model, steady_state):
        # Both auxiliary blocks are built once and stay resident for the whole
        # sweep, because building one costs two solver updates and four name
        # lookups per reaction inside optlang while moving the objective onto
        # it costs one attribute write per reaction. On the 104-organism hCOM
        # community that construction takes around three hours, which the
        # rebuilding version paid on every iteration. A block the objective is
        # not currently on is inert: see :func:`build_push_block` for why.
        #
        # This is worth 1.3-2.5x on models from the textbook up to Recon3D,
        # with the consistent set unchanged in all 19 cases measured. It is
        # not the whole story on very large models: a single build still costs
        # those three hours, which MATLAB does not pay at all because it hands
        # the solver an assembled sparse matrix.
        #
        # Sorted for the same reason the single-shot path sorts: the weights
        # are drawn in one block and zipped positionally, and a set's
        # iteration order varies with PYTHONHASHSEED.
        ordered_ids = sorted(signs)
        reversible = {
            rxn_id
            for rxn_id in ordered_ids
            if is_reversible(model.reactions.get_by_id(rxn_id), signs[rxn_id])
        }
        forward_block = build_push_block(
            model, ordered_ids, signs, tol, reverse=False, prefix="spectra_cc_fwd"
        )
        # The reverse block is built on first use rather than here, because by
        # then the opening forward LP has already pruned everything it can and
        # the block only ever has to cover what is left -- the set of reactions
        # still in play never grows. On a model that is wholly consistent the
        # first pass clears nearly all of it, so deferring is the difference
        # between a block over every reversible reaction and one over a
        # handful.
        reverse_block = None
        model.objective = model.problem.Objective(Zero, direction="max")
        active: Dict["object", float] = {}

        previous_count = None
        while len(rxns_to_check) != previous_count:
            previous_count = len(rxns_to_check)
            n_lps += 2

            batch = [rxn_id for rxn_id in ordered_ids if rxn_id in rxns_to_check]
            model.objective.direction = "max"
            retarget_objective(model, active, draw_weights(forward_block, batch, rng))
            fluxes = _solve(model, "forward consistency LP", False)
            if fluxes is not None:
                rxns_to_check -= carrying_flux(fluxes, detection_cutoff)

            batch = [
                rxn_id
                for rxn_id in ordered_ids
                if rxn_id in rxns_to_check and rxn_id in reversible
            ]
            if reverse_block is None:
                reverse_block = build_push_block(
                    model, batch, signs, tol, reverse=True, prefix="spectra_cc_rev"
                )
            model.objective.direction = "min"
            retarget_objective(model, active, draw_weights(reverse_block, batch, rng))
            fluxes = _solve(model, "reverse LP", False)
            if fluxes is not None:
                rxns_to_check -= carrying_flux(fluxes, detection_cutoff)

            logger.debug(
                "LPs solved: %d - reactions left to check: %d",
                n_lps,
                len(rxns_to_check),
            )

    consistent_ids = set(signs) - rxns_to_check
    logger.info(
        "Final - consistent reactions: %d - inconsistent reactions: %d "
        "[LPs=%d, tol=%.2g, cutoff=%.2g]",
        len(consistent_ids),
        len(rxns_to_check),
        n_lps,
        tol,
        detection_cutoff,
    )
    return consistent_ids, n_lps


def spectra_cc(
    model: "Model",
    tol: float = 1e-4,
    consistency_type: str = STOICHIOMETRY,
    detection_cutoff: Optional[float] = None,
    seed: Optional[int] = None,
    feasibility_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    optimality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
) -> "Model":
    r"""Check the consistency of a metabolic network using SPECTRA.

    Removes the blocked reactions of a model by alternating two LPs, each
    driving as many of the still-unexplained reactions as possible towards the
    flux threshold at once, until that set stops shrinking. Unlike FASTCC no
    reaction flipping is needed: the forward LP's auxiliary variables are
    unbounded from below, so reactions that can only carry negative flux are
    picked up by the reverse LP instead.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on.
    tol : float, optional
        The flux threshold, i.e. the magnitude the LPs drive the reaction
        fluxes to (default 1e-4).
    consistency_type : {"stoichiometry", "topology"}, optional
        Whether to assume a steady state (``S v = 0``) or an accumulation
        condition (``S v >= 0``) (default "stoichiometry"). Topology mode
        treats a reaction as consistent if the metabolites it touches may
        accumulate, which is a weaker requirement, so it never reports fewer
        consistent reactions than stoichiometry mode.
    detection_cutoff : float, optional
        The absolute flux at which a reaction counts as carrying flux
        (default ``0.99 * tol``).
    seed : int, optional
        A seed for the random objective coefficients (default None).
    feasibility_tolerance : float, optional
        How far the mass balance may be violated and still count as
        satisfied. A reaction counts as consistent when an LP drives it to
        ``tol``, so a feasibility tolerance near ``tol`` lets a reaction pass
        on residual rather than on real flux (default 1e-9, the lowest most
        solvers accept). Pass None to leave the solver's own in place.
    optimality_tolerance : float, optional
        How far a reduced cost may be from optimal and still count as optimal
        (default 1e-9). None leaves the solver's own in place.

    Returns
    -------
    cobra.Model
        The consistent model.

    See Also
    --------
    spectra_ccme : Check consistency and extract a model in one pass.

    """
    consistent_ids, _ = consistent_reaction_ids(
        model,
        tol,
        consistency_type,
        detection_cutoff,
        seed,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
    )

    return subset_model(model, consistent_ids)


def blocked_reaction_ids(
    model: "Model",
    tol: float = 1e-4,
    consistency_type: str = STOICHIOMETRY,
    detection_cutoff: Optional[float] = None,
    seed: Optional[int] = None,
    feasibility_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
    optimality_tolerance: Optional[float] = DEFAULT_SOLVER_TOLERANCE,
) -> List[str]:
    """Find the blocked reactions of a model.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on.
    tol : float, optional
        The flux threshold (default 1e-4).
    consistency_type : {"stoichiometry", "topology"}, optional
        Whether to assume a steady state or an accumulation condition
        (default "stoichiometry").
    detection_cutoff : float, optional
        The absolute flux at which a reaction counts as carrying flux
        (default ``0.99 * tol``).
    seed : int, optional
        A seed for the random objective coefficients (default None).
    feasibility_tolerance : float, optional
        How far the mass balance may be violated and still count as satisfied
        (default 1e-9, the lowest most solvers accept). None leaves the
        solver's own in place.
    optimality_tolerance : float, optional
        How far a reduced cost may be from optimal and still count as optimal
        (default 1e-9). None leaves the solver's own in place.

    Returns
    -------
    list of str
        The identifiers of the reactions that cannot carry flux, in the order
        they appear in the model.

    """
    consistent_ids, _ = consistent_reaction_ids(
        model,
        tol,
        consistency_type,
        detection_cutoff,
        seed,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
    )
    return [rxn.id for rxn in model.reactions if rxn.id not in consistent_ids]


def _directions_from_flux(
    fluxes: Dict[str, float], core_ids: Set[str], signs: Dict[str, float]
) -> Dict[str, int]:
    """Read the oriented direction of each core reaction off a flux vector.

    Parameters
    ----------
    fluxes : dict of {str: float}
        The flux through each reaction, as accumulated by the extraction loop.
    core_ids : set of str
        The reactions whose direction matters; everything else gets 0.
    signs : dict of {str: float}
        The orientation of each reaction.

    Returns
    -------
    dict of {str: int}
        1 where the oriented flux is positive, -1 where it is negative, and 0
        for the non-core reactions and any core reaction left at zero flux.

    """
    directions = {}
    for rxn_id, flux in fluxes.items():
        if rxn_id not in core_ids:
            directions[rxn_id] = 0
            continue
        oriented = signs[rxn_id] * flux
        directions[rxn_id] = 1 if oriented > 0 else (-1 if oriented < 0 else 0)
    return directions
