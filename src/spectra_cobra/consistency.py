"""The SPECTRA consistency check."""

from logging import getLogger
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Tuple

import numpy as np

from ._copy import subset_model
from ._lp import carrying_flux, forward_cc, reverse
from ._orientation import (
    STOICHIOMETRY,
    reaction_signs,
    relaxed_mass_balance,
    validate_consistency_type,
)

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)


def consistent_reaction_ids(
    model: "Model",
    tol: float = 1e-4,
    consistency_type: str = STOICHIOMETRY,
    detection_cutoff: Optional[float] = None,
    seed: Optional[int] = None,
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

    signs = reaction_signs(model)
    # Reactions that have not been shown to carry flux yet.
    rxns_to_check = set(signs)
    n_lps = 0

    with model, relaxed_mass_balance(model, steady_state):
        previous_count = None
        while len(rxns_to_check) != previous_count:
            previous_count = len(rxns_to_check)
            n_lps += 2

            with model:
                fluxes = forward_cc(model, rxns_to_check, signs, tol, rng)
            if fluxes is not None:
                rxns_to_check -= carrying_flux(fluxes, detection_cutoff)

            with model:
                fluxes = reverse(model, rxns_to_check, signs, tol, rng)
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

    Returns
    -------
    cobra.Model
        The consistent model.

    See Also
    --------
    spectra_ccme : Check consistency and extract a model in one pass.

    """
    consistent_ids, _ = consistent_reaction_ids(
        model, tol, consistency_type, detection_cutoff, seed
    )

    return subset_model(model, consistent_ids)


def blocked_reaction_ids(
    model: "Model",
    tol: float = 1e-4,
    consistency_type: str = STOICHIOMETRY,
    detection_cutoff: Optional[float] = None,
    seed: Optional[int] = None,
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

    Returns
    -------
    list of str
        The identifiers of the reactions that cannot carry flux, in the order
        they appear in the model.

    """
    consistent_ids, _ = consistent_reaction_ids(
        model, tol, consistency_type, detection_cutoff, seed
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
