"""Finding the smallest community that still does the job.

A large community contains many organisms that contribute nothing to the
function you care about. A *minimal microbiome* is the smallest subset that
still performs it, at a stated fraction of the whole community's growth and
of its production of the metabolites in question.

That is a selection problem over organisms, and it maps onto the extraction
machinery directly. Give each organism's biomass reaction a binary through
``indicator_reactions`` and weight it 1, weight everything else 0, and
minimising the weighted count becomes minimising the number of organisms
present. The biomass coupling built by
:func:`~spectra_cobra.build_community_model` does the rest: an organism
whose binary is off cannot grow, and so carries no flux anywhere.

The growth and production requirements are constraints on the community
model rather than part of the objective, measured first on the full
community and then imposed as a fraction of what it achieved.

    Raghu, A. K., Palanikumar, I., and Raman, K. (2024). Designing
    function-specific minimal microbiomes from large microbial communities.
    *npj Systems Biology and Applications*, 10, 46.
    https://doi.org/10.1038/s41540-024-00373-1
"""

from dataclasses import dataclass
from logging import getLogger
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Tuple

from .community import CommunityModel
from .exceptions import SpectraError
from .formulations import min_net_milp

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)

#: The fraction of the full community's growth a minimal one must keep.
DEFAULT_GROWTH_FRACTION = 0.8

#: The fraction of the full community's production it must keep.
DEFAULT_PRODUCT_FRACTION = 0.8

#: The fraction of each organism's own growth to hold it to while measuring
#: how much the community can produce. Sub-optimal growth is the more
#: realistic setting and gives higher production.
DEFAULT_GROWTH_OPTIMUM_FRACTION = 0.99


@dataclass(frozen=True)
class MinimalMicrobiome:
    """The smallest community found, and what it achieves.

    Parameters
    ----------
    community : CommunityModel
        The minimal community, with the absent organisms removed.
    present : tuple of str
        The organisms retained.
    absent : tuple of str
        The organisms dropped.
    membership : dict of {str: bool}
        Whether each organism of the original community is present.
    growth : float
        The minimal community's total growth.
    production : float
        Its weighted production of the requested metabolites.
    required_growth : float
        The growth it had to reach.
    required_production : float
        The production it had to reach.

    """

    community: CommunityModel
    present: Tuple[str, ...]
    absent: Tuple[str, ...]
    membership: Dict[str, bool]
    growth: float
    production: float
    required_growth: float
    required_production: float

    def summary(self) -> str:
        """Return a one-line account of the result.

        Returns
        -------
        str
            The sizes and what was kept.

        """
        total = len(self.present) + len(self.absent)
        return (
            f"{len(self.present)} of {total} organisms; "
            f"growth {self.growth:.4g} (needed {self.required_growth:.4g}), "
            f"production {self.production:.4g} "
            f"(needed {self.required_production:.4g})"
        )


def _ZERO():
    """Return the solver interface's zero expression.

    Returns
    -------
    optlang expression
        Zero, for building an objective whose terms are set afterwards.

    """
    from optlang.symbolics import Zero

    return Zero


def _linear_coefficients(model: "Model", weights: Dict[str, float]) -> Dict:
    """Return the coefficient of each solver variable in a weighted flux sum.

    Parameters
    ----------
    model : cobra.Model
        The community model.
    weights : dict of {str: float}
        The weight of each reaction.

    Returns
    -------
    dict
        The coefficient of each forward and reverse variable.

    Notes
    -----
    Built as coefficients rather than as a summed expression on purpose.
    Scaling and adding ``flux_expression`` objects produces something the
    solver interface declines to recognise as linear, and it rejects the
    constraint outright.

    """
    coefficients = {}
    for rxn_id, weight in weights.items():
        reaction = model.reactions.get_by_id(rxn_id)
        coefficients[reaction.forward_variable] = float(weight)
        coefficients[reaction.reverse_variable] = -float(weight)
    return coefficients


def _require(model: "Model", weights: Dict[str, float], name: str, lower: float):
    """Add a constraint making a weighted flux sum reach a value.

    Parameters
    ----------
    model : cobra.Model
        The community model.
    weights : dict of {str: float}
        The weight of each reaction.
    name : str
        The constraint's name.
    lower : float
        The value the sum must reach.

    """
    from optlang.symbolics import Zero

    constraint = model.problem.Constraint(Zero, name=name, lb=lower)
    model.add_cons_vars([constraint])
    model.solver.update()
    constraint.set_linear_coefficients(_linear_coefficients(model, weights))


def _weighted_flux(model: "Model", weights: Dict[str, float]) -> float:
    """Return the current value of a weighted flux sum.

    Parameters
    ----------
    model : cobra.Model
        The solved model.
    weights : dict of {str: float}
        The weight of each reaction.

    Returns
    -------
    float
        The weighted sum of the fluxes.

    """
    return sum(
        float(weight) * model.reactions.get_by_id(rxn_id).flux
        for rxn_id, weight in weights.items()
    )


def minimal_microbiome(
    community: CommunityModel,
    products: Optional[Iterable[str]] = None,
    product_weights: Optional[Dict[str, float]] = None,
    growth_fraction: float = DEFAULT_GROWTH_FRACTION,
    product_fraction: float = DEFAULT_PRODUCT_FRACTION,
    growth_optimum_fraction: float = DEFAULT_GROWTH_OPTIMUM_FRACTION,
    required: Optional[Iterable[str]] = None,
    tol: float = 1e-4,
    time_limit: Optional[float] = 300.0,
    seed: Optional[int] = None,
) -> MinimalMicrobiome:
    """Find the smallest subset of a community that still does the job.

    Parameters
    ----------
    community : CommunityModel
        The community to reduce, as built by
        :func:`~spectra_cobra.build_community_model`. Its biomass coupling
        is what makes an absent organism truly absent, so a community built
        with ``couple=False`` is rejected in spirit if not in code: the
        result will keep organisms that contribute only phantom flux.
    products : iterable of str, optional
        The community exchange reactions whose production must be kept.
        With none given only the growth requirement applies.
    product_weights : dict of {str: float}, optional
        A weight per product (default 1 each), for when one product
        matters more than another.
    growth_fraction : float, optional
        The fraction of the full community's growth to keep (default 0.8).
    product_fraction : float, optional
        The fraction of its production to keep (default 0.8).
    growth_optimum_fraction : float, optional
        While measuring how much the full community can produce, hold each
        organism to this fraction of its own growth (default 0.99).
    required : iterable of str, optional
        Organisms that must be kept whatever the objective says.
    tol : float, optional
        The extraction's flux threshold (default 1e-4).
    time_limit : float, optional
        Seconds to spend on the mixed-integer solve (default 300).
    seed : int, optional
        A seed for the extraction's randomised coefficients.

    Returns
    -------
    MinimalMicrobiome
        The minimal community and what it achieves.

    Raises
    ------
    SpectraError
        If the community cannot grow at all, or a named product or
        organism is not part of it.

    Notes
    -----
    The requirements are measured before they are imposed. The full
    community is solved for its growth, then for how much it can produce
    while each organism keeps `growth_optimum_fraction` of its own growth,
    and the minimal community is held to a fraction of each. Asking for a
    fraction of something unmeasured would make the result depend on the
    units of the biomass reaction.

    """
    import numpy as np

    model = community.model
    biomass_ids = [community.biomass_reactions[o] for o in community.organisms]

    required_set = set(required or ())
    unknown = required_set - set(community.organisms)
    if unknown:
        raise SpectraError(f"Not organisms of this community: {sorted(unknown)}.")

    weights: Dict[str, float] = {}
    for rxn_id in products or ():
        if rxn_id not in model.reactions:
            raise SpectraError(f"No reaction {rxn_id!r} in the community model.")
        weights[rxn_id] = float((product_weights or {}).get(rxn_id, 1.0))

    # 1. What the full community manages.
    with model:
        model.objective = {model.reactions.get_by_id(r): 1.0 for r in biomass_ids}
        full_growth = model.slim_optimize()
        if full_growth is None or np.isnan(full_growth) or full_growth <= 0:
            raise SpectraError(
                "The full community cannot grow, so there is no minimal "
                "subset of it that can. Check the medium."
            )
        per_organism = {
            o: model.reactions.get_by_id(community.biomass_reactions[o]).flux
            for o in community.organisms
        }

    # 2. How much it can produce while still growing.
    full_production = 0.0
    if weights:
        with model:
            for organism, value in per_organism.items():
                reaction = model.reactions.get_by_id(
                    community.biomass_reactions[organism]
                )
                reaction.lower_bound = max(
                    reaction.lower_bound, growth_optimum_fraction * value
                )
            from optlang.symbolics import Zero

            model.objective = model.problem.Objective(Zero, direction="max")
            model.objective.set_linear_coefficients(
                _linear_coefficients(model, weights)
            )
            value = model.slim_optimize()
            full_production = 0.0 if value is None or np.isnan(value) else value

    required_growth = growth_fraction * full_growth
    required_production = product_fraction * full_production
    logger.info(
        "full community: growth %.4g, production %.4g; requiring %.4g and %.4g",
        full_growth,
        full_production,
        required_growth,
        required_production,
    )

    # 3. The smallest membership that still meets both.
    with model:
        _require(
            model,
            {r: 1.0 for r in biomass_ids},
            "spectra_community_growth",
            required_growth,
        )
        if weights:
            _require(
                model, weights, "spectra_community_production", required_production
            )
        for organism in required_set:
            reaction = model.reactions.get_by_id(community.biomass_reactions[organism])
            reaction.lower_bound = max(reaction.lower_bound, tol)

        # One binary per organism, on its biomass reaction, and nothing else.
        selection_weights = {
            rxn.id: (1.0 if rxn.id in set(biomass_ids) else 0.0)
            for rxn in model.reactions
        }
        indicators = [
            community.biomass_reactions[o]
            for o in community.organisms
            if o not in required_set
        ]
        # min_net_milp directly rather than through spectra_me: with no core
        # reactions the direction phase has nothing to settle, so the two are
        # the same solve, and this one hands back the indicator values. Those
        # are the membership vector, which is the whole answer here.
        solution = min_net_milp(
            model,
            {},
            selection_weights,
            tol,
            True,
            time_limit,
            None,
            indicators,
        )
        chosen = {
            organism
            for organism in community.organisms
            if community.biomass_reactions[organism] in solution.selected
        }
        present: List[str] = [
            organism
            for organism in community.organisms
            if organism in required_set or organism in chosen
        ]

    membership = {o: (o in set(present)) for o in community.organisms}
    absent = tuple(o for o in community.organisms if not membership[o])
    minimal = community.subset(present)

    # 4. What the minimal community actually achieves.
    reduced = minimal.model
    with reduced:
        reduced.objective = {
            reduced.reactions.get_by_id(minimal.biomass_reactions[o]): 1.0
            for o in minimal.organisms
        }
        growth = reduced.slim_optimize() or 0.0
    production = 0.0
    if weights:
        # Measured the way the requirement was: the most it can produce while
        # still growing. Reading production off the growth optimum instead
        # would report whatever the solver happened to do with the slack,
        # which is usually nothing.
        with reduced:
            kept = {r: w for r, w in weights.items() if r in reduced.reactions}
            _require(
                reduced,
                {minimal.biomass_reactions[o]: 1.0 for o in minimal.organisms},
                "spectra_check_growth",
                required_growth,
            )
            reduced.objective = reduced.problem.Objective(_ZERO(), direction="max")
            reduced.objective.set_linear_coefficients(
                _linear_coefficients(reduced, kept)
            )
            value = reduced.slim_optimize()
            production = 0.0 if value is None or np.isnan(value) else float(value)

    return MinimalMicrobiome(
        community=minimal,
        present=tuple(present),
        absent=absent,
        membership=membership,
        growth=float(growth),
        production=float(production),
        required_growth=float(required_growth),
        required_production=float(required_production),
    )
