"""Joining organism models into a community model.

A community model puts every organism's network into its own compartment and
lets them trade through a shared pool. Each organism keeps its own
extracellular space; a transport reaction moves each exchangeable metabolite
between that space and the shared pool, and a single community exchange
connects the pool to the environment. One organism's secretion is therefore
available to another, which is what makes cross-feeding possible.

Two things need care. An organism that is "absent" has to be absent in full,
not merely not growing, so every reaction is coupled to its own biomass: no
biomass flux, no flux anywhere in that organism. And the environment must be
reachable only through the shared pool, so the organisms' own exchange
reactions are replaced rather than kept alongside.

The compartmentalised layout and the biomass coupling follow the joint FBA
formulation used for minimal microbiomes:

    Raghu, A. K., Palanikumar, I., and Raman, K. (2024). Designing
    function-specific minimal microbiomes from large microbial communities.
    *npj Systems Biology and Applications*, 10, 46.
    https://doi.org/10.1038/s41540-024-00373-1
"""

from dataclasses import dataclass, field
from logging import getLogger
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Sequence, Tuple

from .exceptions import SpectraError

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)

#: Separator between a reaction or metabolite's own identifier and the
#: organism it belongs to.
ORGANISM_SEPARATOR = "__"

#: Suffix marking the shared pool a community trades through.
SHARED_SUFFIX = "u"

#: How much flux a reaction may carry per unit of its organism's growth.
#: Large enough not to bind on a growing organism, small enough that an
#: absent one carries nothing.
COUPLING_FACTOR = 1000.0

#: The flux a coupled reaction may carry even at zero growth, which keeps
#: the constraint from forcing an exactly-zero solution the solver may
#: struggle to find.
COUPLING_THRESHOLD = 0.01


@dataclass(frozen=True)
class CommunityModel:
    """A community of organism models and the bookkeeping to read it.

    Parameters
    ----------
    model : cobra.Model
        The joined model.
    organisms : tuple of str
        The organism identifiers, in the order given.
    biomass_reactions : dict of {str: str}
        The biomass reaction of each organism, keyed by organism.
    reactions_of : dict of {str: tuple of str}
        The reactions belonging to each organism, keyed by organism.
    community_exchanges : tuple of str
        The exchange reactions connecting the shared pool to the
        environment.

    """

    model: "Model"
    organisms: Tuple[str, ...]
    biomass_reactions: Dict[str, str]
    reactions_of: Dict[str, Tuple[str, ...]] = field(repr=False, default_factory=dict)
    community_exchanges: Tuple[str, ...] = ()

    def organism_of(self, reaction_id: str) -> Optional[str]:
        """Return the organism a reaction belongs to, if any.

        Parameters
        ----------
        reaction_id : str
            The reaction to look up.

        Returns
        -------
        str or None
            The organism, or None for a community exchange.

        """
        if ORGANISM_SEPARATOR not in reaction_id:
            return None
        candidate = reaction_id.rsplit(ORGANISM_SEPARATOR, 1)[-1]
        return candidate if candidate in self.biomass_reactions else None

    def subset(self, keep: Iterable[str]) -> "CommunityModel":
        """Return the community restricted to the given organisms.

        Parameters
        ----------
        keep : iterable of str
            The organisms to retain.

        Returns
        -------
        CommunityModel
            A new community holding only those organisms. The shared pool
            and its exchanges are kept, so the survivors still trade and
            still face the same environment.

        Raises
        ------
        SpectraError
            If an organism is named that the community does not have.

        """
        wanted = list(dict.fromkeys(keep))
        unknown = set(wanted) - set(self.organisms)
        if unknown:
            raise SpectraError(
                f"Not organisms of this community: {sorted(unknown)[:5]}."
            )
        dropped = [org for org in self.organisms if org not in wanted]
        model = self.model.copy()
        model.remove_reactions(
            [r for org in dropped for r in self.reactions_of[org]],
            remove_orphans=True,
        )
        return CommunityModel(
            model=model,
            organisms=tuple(org for org in self.organisms if org in wanted),
            biomass_reactions={
                k: v for k, v in self.biomass_reactions.items() if k in wanted
            },
            reactions_of={k: v for k, v in self.reactions_of.items() if k in wanted},
            community_exchanges=self.community_exchanges,
        )


def _find_biomass(model: "Model") -> str:
    """Return a model's biomass reaction.

    Parameters
    ----------
    model : cobra.Model
        The model to inspect.

    Returns
    -------
    str
        The reaction identifier.

    Raises
    ------
    SpectraError
        If no biomass reaction can be identified.

    """
    objective = [r.id for r in model.reactions if r.objective_coefficient != 0]
    if len(objective) == 1:
        return objective[0]
    named = [r.id for r in model.reactions if r.id.lower().startswith("biomass")]
    if len(named) == 1:
        return named[0]
    raise SpectraError(
        f"Cannot tell which reaction of {model.id!r} is its biomass reaction "
        f"({len(objective)} in the objective, {len(named)} named 'biomass'). "
        f"Pass biomass_reactions explicitly."
    )


def _couple_to_biomass(
    model: "Model",
    reaction_ids: Sequence[str],
    biomass_id: str,
    factor: float,
    threshold: float,
) -> List["object"]:
    """Build the constraints tying an organism's reactions to its growth.

    Parameters
    ----------
    model : cobra.Model
        The community model.
    reaction_ids : sequence of str
        The reactions to couple.
    biomass_id : str
        The organism's biomass reaction.
    factor : float
        How much flux is allowed per unit of growth.
    threshold : float
        The flux allowed at zero growth.

    Returns
    -------
    list
        The constraints to add.

    Notes
    -----
    Two inequalities per reaction, ``v - c*b <= u`` and ``v + c*b >= -u``,
    so a reaction can carry flux in either direction only in proportion to
    its organism's growth. Without them an absent organism could still run
    its metabolism and feed the rest of the community.

    """
    prob = model.problem
    biomass = model.reactions.get_by_id(biomass_id).flux_expression
    constraints = []
    for rxn_id in reaction_ids:
        if rxn_id == biomass_id:
            continue
        flux = model.reactions.get_by_id(rxn_id).flux_expression
        constraints.append(
            prob.Constraint(
                flux - factor * biomass, name=f"spectra_couple_u_{rxn_id}", ub=threshold
            )
        )
        constraints.append(
            prob.Constraint(
                flux + factor * biomass,
                name=f"spectra_couple_l_{rxn_id}",
                lb=-threshold,
            )
        )
    return constraints


def build_community_model(
    models: Sequence["Model"],
    organisms: Optional[Sequence[str]] = None,
    biomass_reactions: Optional[Dict[str, str]] = None,
    shared_compartment: str = SHARED_SUFFIX,
    external_compartment: str = "e",
    couple: bool = True,
    coupling_factor: float = COUPLING_FACTOR,
    coupling_threshold: float = COUPLING_THRESHOLD,
    model_id: str = "community",
) -> CommunityModel:
    """Join organism models into a community that trades through a shared pool.

    Parameters
    ----------
    models : sequence of cobra.Model
        The organism models. They are not modified.
    organisms : sequence of str, optional
        A short identifier for each, in the same order (default: each
        model's own id). These become the suffix on every reaction and
        metabolite, so they must be unique and free of the separator
        ``__``.
    biomass_reactions : dict of {str: str}, optional
        The biomass reaction of each organism, keyed by organism (default:
        the objective, or a reaction named ``biomass*``).
    shared_compartment : str, optional
        The compartment of the shared pool (default ``"u"``).
    external_compartment : str, optional
        The compartment whose metabolites are exchangeable (default
        ``"e"``).
    couple : bool, optional
        Whether to tie each organism's reactions to its own biomass
        (default True). Without this an organism that is not growing can
        still run its metabolism and feed the others.
    coupling_factor : float, optional
        Flux allowed per unit of growth (default 1000).
    coupling_threshold : float, optional
        Flux allowed at zero growth (default 0.01).
    model_id : str, optional
        The identifier of the joined model (default "community").

    Returns
    -------
    CommunityModel
        The joined model and the bookkeeping needed to read it.

    Raises
    ------
    SpectraError
        If the organism names are not usable, or a biomass reaction cannot
        be identified.

    Notes
    -----
    Each organism keeps its own copy of the external compartment, and a
    transport reaction links it to the shared pool, so metabolites move
    organism -> pool -> organism. The organisms' own exchange reactions are
    dropped: leaving them in place would let every organism draw on the
    environment directly, and no amount of community structure would
    constrain anything.

    """
    from cobra import Metabolite, Model, Reaction

    if organisms is None:
        organisms = [model.id for model in models]
    if len(organisms) != len(models):
        raise SpectraError(
            f"Got {len(models)} models but {len(organisms)} organism names."
        )
    if len(set(organisms)) != len(organisms):
        raise SpectraError("Organism names must be unique.")
    bad = [o for o in organisms if ORGANISM_SEPARATOR in o or not o]
    if bad:
        raise SpectraError(
            f"Organism names must be non-empty and must not contain "
            f"{ORGANISM_SEPARATOR!r}: {bad[:5]}."
        )

    resolved = dict(biomass_reactions or {})
    for model, organism in zip(models, organisms):
        if organism not in resolved:
            resolved[organism] = _find_biomass(model)

    community = Model(model_id)
    shared: Dict[str, "Metabolite"] = {}
    reactions_of: Dict[str, Tuple[str, ...]] = {}
    biomass_ids: Dict[str, str] = {}
    transports: List["Reaction"] = []

    def tag(identifier: str, organism: str) -> str:
        return f"{identifier}{ORGANISM_SEPARATOR}{organism}"

    for model, organism in zip(models, organisms):
        # Only the exchanges facing the environment are replaced. A demand
        # or sink inside the cell is also a boundary reaction by cobrapy's
        # reckoning, and a biomass reaction that only consumes looks exactly
        # like one, so dropping every boundary reaction would delete it.
        exchanges = {
            rxn.id
            for rxn in model.boundary
            if any(met.compartment == external_compartment for met in rxn.metabolites)
        }
        metabolites = {
            met.id: Metabolite(
                tag(met.id, organism),
                formula=met.formula,
                name=met.name,
                charge=met.charge,
                compartment=f"{met.compartment}{ORGANISM_SEPARATOR}{organism}",
            )
            for met in model.metabolites
        }
        community.add_metabolites(list(metabolites.values()))

        own: List["Reaction"] = []
        for reaction in model.reactions:
            if reaction.id in exchanges:
                # Replaced by a transport into the shared pool below.
                continue
            copied = Reaction(
                tag(reaction.id, organism),
                name=reaction.name,
                lower_bound=reaction.lower_bound,
                upper_bound=reaction.upper_bound,
            )
            own.append((copied, reaction))
        community.add_reactions([c for c, _ in own])
        for copied, reaction in own:
            copied.add_metabolites(
                {metabolites[m.id]: c for m, c in reaction.metabolites.items()}
            )

        # One transport per exchangeable metabolite, organism <-> pool.
        for met in model.metabolites:
            if met.compartment != external_compartment:
                continue
            if met.id not in shared:
                # Both conventions for naming a compartment are in use:
                # "glc_D_e" and "glc_D[e]". Strip either so the shared pool
                # is named after the metabolite rather than after one
                # organism's spelling of where it sits.
                base = met.id
                for suffix in (
                    f"_{external_compartment}",
                    f"[{external_compartment}]",
                ):
                    if base.endswith(suffix):
                        base = base[: -len(suffix)]
                        break
                pooled = Metabolite(
                    f"{base}_{shared_compartment}",
                    formula=met.formula,
                    name=met.name,
                    compartment=shared_compartment,
                    charge=met.charge,
                )
                shared[met.id] = pooled
                community.add_metabolites([pooled])
            transport = Reaction(
                tag(f"TR_{met.id}", organism),
                name=f"{met.name} transport, {organism}",
                lower_bound=-1000.0,
                upper_bound=1000.0,
            )
            transports.append(
                (transport, organism, metabolites[met.id], shared[met.id])
            )

        reactions_of[organism] = tuple(c.id for c, _ in own)
        biomass_ids[organism] = tag(resolved[organism], organism)

    community.add_reactions([t for t, _, _, _ in transports])
    for transport, organism, local, pooled in transports:
        transport.add_metabolites({local: -1.0, pooled: 1.0})
        reactions_of[organism] = reactions_of[organism] + (transport.id,)

    # The environment, reachable only through the pool.
    community_exchanges = []
    for pooled in shared.values():
        exchange = Reaction(
            f"EX_{pooled.id}",
            name=f"{pooled.name} community exchange",
            lower_bound=0.0,
            upper_bound=1000.0,
        )
        community.add_reactions([exchange])
        exchange.add_metabolites({pooled: -1.0})
        community_exchanges.append(exchange.id)

    missing = [o for o in organisms if biomass_ids[o] not in community.reactions]
    if missing:
        raise SpectraError(
            f"The biomass reaction of {missing[:5]} is not in the community; "
            f"it may have been an exchange reaction, which is not kept."
        )

    community.objective = {
        community.reactions.get_by_id(biomass_ids[o]): 1.0 for o in organisms
    }

    if couple:
        constraints = []
        for organism in organisms:
            constraints.extend(
                _couple_to_biomass(
                    community,
                    reactions_of[organism],
                    biomass_ids[organism],
                    coupling_factor,
                    coupling_threshold,
                )
            )
        community.add_cons_vars(constraints)

    logger.info(
        "community of %d organisms: %d reactions, %d metabolites, "
        "%d shared metabolites",
        len(organisms),
        len(community.reactions),
        len(community.metabolites),
        len(shared),
    )
    return CommunityModel(
        model=community,
        organisms=tuple(organisms),
        biomass_reactions=biomass_ids,
        reactions_of=reactions_of,
        community_exchanges=tuple(community_exchanges),
    )
