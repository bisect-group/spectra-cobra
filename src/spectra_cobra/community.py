"""Joining models into one network whose parts trade with each other.

The parts may be organisms of a microbial community or tissues of a body.
Nothing in the construction distinguishes them, so this module calls them
*units* and names them `organisms` in the signatures for continuity.

Every unit keeps its own copy of the network, tagged so that its reactions
and metabolites cannot be confused with another's, and the units meet in a
**pool**: a compartment several of them can reach. One unit's secretion
arrives in the pool and is available to the rest, which is what makes
cross-feeding possible. A pool named in `environment` also carries exchange
reactions, which is how anything enters or leaves the whole system.

There are two ways to build the pool, and both are in published use:

``"pooled"``
    Each unit keeps a private external compartment, and a transport
    reaction moves each exchangeable metabolite between it and the pool.
    Because the transport belongs to one unit, it can be bounded per unit
    -- which is where measured uptake and secretion rates go. Several pools
    can coexist, so the units need not all meet in the same place.

``"shared"``
    The units share the external compartment outright: one row of the
    stoichiometric matrix per exchangeable metabolite, for the whole
    community. There are no transport reactions and there is only ever one
    pool. For a hundred units across a few thousand metabolites this saves
    a hundred thousand reactions, which is the difference between a model
    that can be built and one that cannot.

The compartmentalised layout, the biomass coupling and its parameters
follow the joint FBA formulation used for minimal microbiomes:

    Raghu, A. K., Palanikumar, I., and Raman, K. (2024). Designing
    function-specific minimal microbiomes from large microbial communities.
    *npj Systems Biology and Applications*, 10, 46.
    https://doi.org/10.1038/s41540-024-00373-1
"""

from dataclasses import dataclass, field
from logging import getLogger
from typing import (
    TYPE_CHECKING,
    Callable,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
)

from ._copy import subset_model
from .exceptions import SpectraError

if TYPE_CHECKING:
    from cobra.core import Metabolite, Model


logger = getLogger(__name__)

#: Separator between a reaction or metabolite's own identifier and the unit
#: it belongs to.
ORGANISM_SEPARATOR = "__"

#: Name of the pool the units trade through, when only one is needed.
SHARED_SUFFIX = "u"

#: Each unit keeps its own external compartment and reaches the pool
#: through transport reactions.
POOLED = "pooled"

#: The units share the external compartment itself, with no transports.
SHARED = "shared"

#: The ways of building a pool.
POOL_MODES = (POOLED, SHARED)

#: How much flux a reaction may carry per unit of its anchor's flux. The
#: default is the published value for a biomass anchor on a microbial
#: community; an anchor that runs at a different scale needs a different
#: factor, and :meth:`CommunityModel.coupling_summary` reports what the
#: chosen one amounts to.
COUPLING_FACTOR = 1000.0

#: The flux a coupled reaction may carry at zero anchor flux, which keeps
#: the constraint from forcing an exactly-zero solution the solver may
#: struggle to find.
COUPLING_THRESHOLD = 0.01


@dataclass(frozen=True)
class CommunityModel:
    """Units joined into one model, and the bookkeeping to read it.

    Parameters
    ----------
    model : cobra.Model
        The joined model.
    organisms : tuple of str
        The unit identifiers, in the order given.
    biomass_reactions : dict of {str: str}
        The anchor reaction of each unit, keyed by unit. Named for the
        usual case, in which the anchor is a biomass reaction.
    reactions_of : dict of {str: tuple of str}
        The reactions belonging to each unit, keyed by unit.
    community_exchanges : tuple of str
        The exchange reactions connecting a pool to the environment.
    mode : str
        How the pool was built, ``"pooled"`` or ``"shared"``.
    pools : dict of {str: tuple of str}
        The units meeting in each pool, keyed by pool.
    pool_metabolites : dict of {str: tuple of str}
        The metabolites of each pool, keyed by pool.
    anchor_sources : dict of {str: str}
        How each unit's anchor was decided: ``"given"``, ``"objective"``
        or ``"name"``.
    anchor_capacity : dict of {str: float}
        The most flux each unit's anchor can carry in the community before
        coupling is imposed, or empty if the check was skipped.
    draft_reactions : dict of {str: tuple of str}
        Where a database was supplied, the reactions that came from the
        unit's own model rather than from the database.
    database_reactions : dict of {str: tuple of str}
        The reactions that came only from the database, and are therefore
        candidates to be gap-filled in.
    coupling_factor : float or None
        The factor used, or None if the units were not coupled.
    coupling_threshold : float or None
        The threshold used, or None if the units were not coupled.

    """

    model: "Model"
    organisms: Tuple[str, ...]
    biomass_reactions: Dict[str, str]
    reactions_of: Dict[str, Tuple[str, ...]] = field(repr=False, default_factory=dict)
    community_exchanges: Tuple[str, ...] = ()
    mode: str = POOLED
    pools: Dict[str, Tuple[str, ...]] = field(default_factory=dict)
    pool_metabolites: Dict[str, Tuple[str, ...]] = field(
        repr=False, default_factory=dict
    )
    anchor_sources: Dict[str, str] = field(default_factory=dict)
    anchor_capacity: Dict[str, float] = field(default_factory=dict)
    draft_reactions: Dict[str, Tuple[str, ...]] = field(
        repr=False, default_factory=dict
    )
    database_reactions: Dict[str, Tuple[str, ...]] = field(
        repr=False, default_factory=dict
    )
    coupling_factor: Optional[float] = None
    coupling_threshold: Optional[float] = None

    @property
    def anchor_reactions(self) -> Dict[str, str]:
        """Return the anchor reaction of each unit.

        Returns
        -------
        dict of {str: str}
            The same mapping as `biomass_reactions`, under the name that
            fits a unit whose anchor is not a biomass reaction.

        """
        return self.biomass_reactions

    def organism_of(self, reaction_id: str) -> Optional[str]:
        """Return the unit a reaction belongs to, if any.

        Parameters
        ----------
        reaction_id : str
            The reaction to look up.

        Returns
        -------
        str or None
            The unit, or None for a community exchange.

        """
        if ORGANISM_SEPARATOR not in reaction_id:
            return None
        candidate = reaction_id.rsplit(ORGANISM_SEPARATOR, 1)[-1]
        return candidate if candidate in self.biomass_reactions else None

    def coupling_summary(self) -> str:
        """Return a table of what each unit is coupled to.

        Returns
        -------
        str
            One line per unit, giving its anchor, how the anchor was
            decided, and -- where it was measured -- the most flux the
            anchor can carry and the cap that imposes on the unit.

        Notes
        -----
        The cap is ``coupling_factor * capacity + coupling_threshold``: the
        most flux any reaction of the unit may carry when its anchor runs
        as fast as it can. A cap far below the unit's own bounds means the
        coupling is throttling it, and the factor wants raising.

        """
        if not self.coupling_factor:
            return "not coupled"
        rows = [f"{'unit':<16} {'anchor':<32} {'how':<10} {'capacity':>10} {'cap':>10}"]
        for organism in self.organisms:
            capacity = self.anchor_capacity.get(organism)
            if capacity is None:
                shown, cap = "-", "-"
            else:
                shown = f"{capacity:.4g}"
                limit = self.coupling_factor * capacity + (
                    self.coupling_threshold or 0.0
                )
                cap = f"{limit:.4g}"
            rows.append(
                f"{organism:<16} {self.biomass_reactions[organism]:<32} "
                f"{self.anchor_sources.get(organism, '-'):<10} {shown:>10} {cap:>10}"
            )
        return "\n".join(rows)

    def subset(self, keep: Iterable[str]) -> "CommunityModel":
        """Return the community restricted to the given units.

        Parameters
        ----------
        keep : iterable of str
            The units to retain.

        Returns
        -------
        CommunityModel
            A new community holding only those units. The pools and their
            exchanges are kept, so the survivors still trade and still
            face the same environment.

        Raises
        ------
        SpectraError
            If a unit is named that the community does not have.

        """
        wanted = list(dict.fromkeys(keep))
        unknown = set(wanted) - set(self.organisms)
        if unknown:
            raise SpectraError(
                f"Not organisms of this community: {sorted(unknown)[:5]}."
            )
        dropped = [org for org in self.organisms if org not in wanted]
        leaving = {r for org in dropped for r in self.reactions_of[org]}
        model = subset_model(
            self.model, {r.id for r in self.model.reactions} - leaving
        )
        survivors = tuple(org for org in self.organisms if org in wanted)

        def only(mapping):
            return {k: v for k, v in mapping.items() if k in set(wanted)}

        return CommunityModel(
            model=model,
            organisms=survivors,
            biomass_reactions=only(self.biomass_reactions),
            reactions_of=only(self.reactions_of),
            community_exchanges=self.community_exchanges,
            mode=self.mode,
            pools={
                pool: tuple(o for o in members if o in set(wanted))
                for pool, members in self.pools.items()
            },
            pool_metabolites=self.pool_metabolites,
            anchor_sources=only(self.anchor_sources),
            anchor_capacity=only(self.anchor_capacity),
            draft_reactions=only(self.draft_reactions),
            database_reactions=only(self.database_reactions),
            coupling_factor=self.coupling_factor,
            coupling_threshold=self.coupling_threshold,
        )

    def with_model(self, model: "Model") -> "CommunityModel":
        """Return the same community over a reduced model.

        Parameters
        ----------
        model : cobra.Model
            A model holding some subset of this community's reactions,
            as a consistency check or an extraction leaves behind.

        Returns
        -------
        CommunityModel
            The bookkeeping narrowed to what survived, so that a unit's
            reactions, its database candidates and the pool exchanges all
            still name reactions that exist.

        """
        present = {rxn.id for rxn in model.reactions}

        def kept(mapping):
            return {
                unit: tuple(r for r in ids if r in present)
                for unit, ids in mapping.items()
            }

        return CommunityModel(
            model=model,
            organisms=self.organisms,
            biomass_reactions=dict(self.biomass_reactions),
            reactions_of=kept(self.reactions_of),
            community_exchanges=tuple(
                r for r in self.community_exchanges if r in present
            ),
            mode=self.mode,
            pools=dict(self.pools),
            pool_metabolites={
                pool: tuple(
                    met for met in mets if met in {m.id for m in model.metabolites}
                )
                for pool, mets in self.pool_metabolites.items()
            },
            anchor_sources=dict(self.anchor_sources),
            anchor_capacity=dict(self.anchor_capacity),
            draft_reactions=kept(self.draft_reactions),
            database_reactions=kept(self.database_reactions),
            coupling_factor=self.coupling_factor,
            coupling_threshold=self.coupling_threshold,
        )

    def decompose(self) -> Dict[str, "Model"]:
        """Split the community back into one model per unit.

        Returns
        -------
        dict of {str: cobra.Model}
            Each unit's own network, with the tags stripped and an
            exchange reaction wherever the unit met a pool.

        Notes
        -----
        What was a transport into the pool becomes an exchange: dropping
        the pool side of ``x[e] -> x_u`` leaves ``x[e] ->``, which is an
        exchange carrying the transport's own bounds. Where the units
        shared the external compartment there are no transports, so an
        exchange is added for every external metabolite the unit touches.

        """
        from cobra import Model

        pooled_ids = {
            met_id for mets in self.pool_metabolites.values() for met_id in mets
        }
        models: Dict[str, "Model"] = {}
        for organism in self.organisms:
            tagged = set(self.reactions_of[organism])
            extracted = subset_model(self.model, tagged)
            suffix = f"{ORGANISM_SEPARATOR}{organism}"
            # The pool is shared, so it cannot travel with any one unit.
            # Dropping it turns each transport into an exchange.
            dangling = [m for m in extracted.metabolites if m.id in pooled_ids]
            if self.mode == SHARED:
                for met in dangling:
                    if met.reactions:
                        extracted.add_boundary(met, type="exchange")
            else:
                # Reversed for the same reason as elsewhere: deleting
                # from the front of optlang's container re-indexes
                # everything after it.
                extracted.remove_metabolites(
                    list(reversed(dangling)), destructive=False
                )

            own = Model(organism)
            own.solver = self.model.problem
            own.add_reactions([r.copy() for r in extracted.reactions])
            for reaction in own.reactions:
                if reaction.id.endswith(suffix):
                    reaction.id = reaction.id[: -len(suffix)]
            for metabolite in own.metabolites:
                if metabolite.id.endswith(suffix):
                    metabolite.id = metabolite.id[: -len(suffix)]
                if metabolite.compartment and metabolite.compartment.endswith(suffix):
                    metabolite.compartment = metabolite.compartment[: -len(suffix)]
            own.repair()
            anchor = self.biomass_reactions[organism]
            if anchor.endswith(suffix):
                anchor = anchor[: -len(suffix)]
            if anchor in own.reactions:
                own.objective = own.reactions.get_by_id(anchor)
            models[organism] = own
        return models


def _strip_compartment(met_id: str, compartment: str) -> str:
    """Return a metabolite identifier without its compartment suffix.

    Parameters
    ----------
    met_id : str
        The identifier.
    compartment : str
        The compartment to strip.

    Returns
    -------
    str
        The identifier with ``_c`` or ``[c]`` removed, if present.

    Notes
    -----
    Both conventions are in use -- ``glc_D_e`` and ``glc_D[e]`` -- and
    getting this wrong does not fail loudly. Each unit would simply trade
    through its own private pool, named after its own spelling, and the
    community would never cross-feed.

    """
    for suffix in (f"_{compartment}", f"[{compartment}]"):
        if met_id.endswith(suffix):
            return met_id[: -len(suffix)]
    return met_id


def _find_anchor(model: "Model", organism: str) -> Tuple[str, str]:
    """Return the reaction a unit's coupling should gate on.

    Parameters
    ----------
    model : cobra.Model
        The unit's model.
    organism : str
        Its name, for the error message.

    Returns
    -------
    tuple of (str, str)
        The reaction identifier, and how it was decided.

    Raises
    ------
    SpectraError
        If no single reaction stands out.

    Notes
    -----
    Guessing is deliberately avoided. Coupling every reaction of a unit to
    the wrong one would throttle or free the whole unit without any error,
    and a human model's anchor may be growth, maintenance, or a demand
    reaction with no telling name.

    """
    objective = [r.id for r in model.reactions if r.objective_coefficient != 0]
    if len(objective) == 1:
        return objective[0], "objective"
    named = [r.id for r in model.reactions if r.id.lower().startswith("biomass")]
    if len(named) == 1:
        return named[0], "name"
    candidates = sorted(
        r.id
        for r in model.reactions
        if r.id.lower().startswith("biomass")
        or "maintenance" in r.id.lower()
        or r.id.upper() in {"ATPM", "DM_ATP_C_", "DM_ATP_C"}
    )
    raise SpectraError(
        f"Cannot tell which reaction of {organism!r} its parts should be "
        f"coupled to ({len(objective)} in the objective, {len(named)} named "
        f"'biomass'). Pass biomass_reactions explicitly"
        + (f"; candidates look like {candidates[:5]}." if candidates else ".")
    )


def _couple_to_anchor(
    model: "Model",
    reaction_ids: Sequence[str],
    anchor_id: str,
    factor: float,
    threshold: float,
) -> List[object]:
    """Build the constraints tying a unit's reactions to its anchor.

    Parameters
    ----------
    model : cobra.Model
        The community model.
    reaction_ids : sequence of str
        The reactions to couple.
    anchor_id : str
        The unit's anchor reaction.
    factor : float
        How much flux is allowed per unit of anchor flux.
    threshold : float
        The flux allowed at zero anchor flux.

    Returns
    -------
    list
        The constraints to add.

    Notes
    -----
    Two inequalities per reaction, ``v - c*a <= u`` and ``v + c*a >= -u``,
    so a reaction can carry flux in either direction only in proportion to
    its unit's anchor. Without them an absent unit could still run its
    metabolism and feed the rest of the community.

    """
    prob = model.problem
    anchor = model.reactions.get_by_id(anchor_id).flux_expression
    constraints = []
    for rxn_id in reaction_ids:
        if rxn_id == anchor_id:
            continue
        flux = model.reactions.get_by_id(rxn_id).flux_expression
        constraints.append(
            prob.Constraint(
                flux - factor * anchor, name=f"spectra_couple_u_{rxn_id}", ub=threshold
            )
        )
        constraints.append(
            prob.Constraint(
                flux + factor * anchor,
                name=f"spectra_couple_l_{rxn_id}",
                lb=-threshold,
            )
        )
    return constraints


def _resolve_pools(
    organisms: Sequence[str],
    pools: Optional[Mapping[str, Iterable[str]]],
    mode: str,
    shared_compartment: str,
) -> Dict[str, Tuple[str, ...]]:
    """Return the units meeting in each pool.

    Parameters
    ----------
    organisms : sequence of str
        Every unit.
    pools : mapping, optional
        The requested pools, or None for one pool holding everything.
    mode : str
        How the pool is built.
    shared_compartment : str
        The name to give the default pool.

    Returns
    -------
    dict of {str: tuple of str}
        The members of each pool.

    Raises
    ------
    SpectraError
        If the pools name unknown units, leave a unit out, or ask for
        several pools in a mode that can only have one.

    """
    if pools is None:
        return {shared_compartment: tuple(organisms)}
    resolved = {pool: tuple(dict.fromkeys(m)) for pool, m in pools.items()}
    if not resolved:
        raise SpectraError("pools must name at least one pool.")
    unknown = {o for members in resolved.values() for o in members} - set(organisms)
    if unknown:
        raise SpectraError(
            f"pools name units that do not exist: {sorted(unknown)[:5]}."
        )
    homeless = [o for o in organisms if not any(o in m for m in resolved.values())]
    if homeless:
        raise SpectraError(
            f"Every unit must be in at least one pool, but {homeless[:5]} is "
            f"in none. A unit in no pool cannot exchange anything and cannot "
            f"grow."
        )
    if mode == SHARED and len(resolved) > 1:
        raise SpectraError(
            f"mode={SHARED!r} shares the external compartment itself, so "
            f"there can only be one pool; got {len(resolved)}. Use "
            f"mode={POOLED!r} for units that must meet in different places."
        )
    clash = set(resolved) & set(organisms)
    if clash:
        raise SpectraError(f"A pool cannot be named after a unit: {sorted(clash)}.")
    return resolved


def _resolve_environment(
    environment: Optional[Iterable[str]], pools: Mapping[str, Sequence[str]]
) -> Tuple[str, ...]:
    """Return the pools that reach the environment.

    Parameters
    ----------
    environment : iterable of str, optional
        The requested pools, or None.
    pools : mapping
        Every pool.

    Returns
    -------
    tuple of str
        The pools to give exchange reactions.

    Raises
    ------
    SpectraError
        If an unknown pool is named, or none is named and there is more
        than one to choose from.

    Notes
    -----
    Defaulting to every pool would be wrong. A pool that exists only to
    link two units is an interface, and giving it exchanges would let
    those two dump metabolites straight into the environment, bypassing
    whichever pool was supposed to be the only way out.

    """
    if environment is None:
        if len(pools) == 1:
            return tuple(pools)
        raise SpectraError(
            f"With {len(pools)} pools, environment must say which of them "
            f"reach the environment: {sorted(pools)}. A pool that only links "
            f"units together must stay closed, or those units could bypass "
            f"the rest of the system."
        )
    wanted = tuple(dict.fromkeys(environment))
    unknown = set(wanted) - set(pools)
    if unknown:
        raise SpectraError(f"environment names unknown pools: {sorted(unknown)}.")
    return wanted


def _merge_with_database(
    draft: "Model", database: "Model", organism: str
) -> Tuple["Model", List[str], List[str]]:
    """Return a unit holding a whole database, bounded by its own model.

    Parameters
    ----------
    draft : cobra.Model
        The unit's own reconstruction.
    database : cobra.Model
        The reactions it may draw on.
    organism : str
        Its name, for the log line.

    Returns
    -------
    tuple of (cobra.Model, list of str, list of str)
        The merged model, the reactions the draft already had, and the
        reactions only the database has.

    Notes
    -----
    The draft contributes its bounds rather than its reactions: where an
    identifier appears in both, the draft's bounds win, because they are
    what is known about this organism. Reactions the draft has and the
    database does not are kept, which a draft carved out of the database
    will rarely have beyond its biomass reaction.

    """
    from cobra import Metabolite, Model, Reaction

    database_ids = {r.id for r in database.reactions}
    draft_ids = {r.id for r in draft.reactions}
    bounds = {r.id: r.bounds for r in draft.reactions}

    # Built from fresh objects rather than by copying either model.
    # cobrapy's copy deep-copies the solver, which for a universal
    # reconstruction is slow and, on some solvers, simply fails -- and
    # neither input is modified here, so there is nothing to protect.
    merged = Model(organism)
    metabolites: Dict[str, "Metabolite"] = {}
    for source in (database, draft):
        for met in source.metabolites:
            if met.id not in metabolites:
                metabolites[met.id] = Metabolite(
                    met.id,
                    formula=met.formula,
                    name=met.name,
                    charge=met.charge,
                    compartment=met.compartment,
                )
    merged.add_metabolites(list(metabolites.values()))

    extra = [r for r in draft.reactions if r.id not in database_ids]
    copied = []
    for reaction in list(database.reactions) + extra:
        lower, upper = bounds.get(reaction.id, reaction.bounds)
        fresh = Reaction(
            reaction.id,
            name=reaction.name,
            subsystem=reaction.subsystem,
            lower_bound=lower,
            upper_bound=upper,
        )
        copied.append((fresh, reaction))
    merged.add_reactions([f for f, _ in copied])
    for fresh, reaction in copied:
        fresh.add_metabolites(
            {metabolites[m.id]: c for m, c in reaction.metabolites.items()}
        )

    objective = [r.id for r in draft.reactions if r.objective_coefficient != 0]
    if len(objective) == 1 and objective[0] in merged.reactions:
        merged.objective = merged.reactions.get_by_id(objective[0])
    logger.info(
        "%s: %d reactions from the database, %d from the draft (%d of them "
        "the database does not have)",
        organism,
        len(database_ids),
        len(draft_ids),
        len(extra),
    )
    return merged, sorted(draft_ids), sorted(database_ids - draft_ids)


def _external_compartments(
    units: Sequence["Model"],
    organisms: Sequence[str],
    given: Optional[Union[str, Mapping[str, str]]],
) -> Dict[str, str]:
    """Return the external compartment of each unit.

    Parameters
    ----------
    units : sequence of cobra.Model
        The models.
    organisms : sequence of str
        Their names, in the same order.
    given : str or mapping, optional
        One compartment for all of them, or one per unit, or None to work
        it out from each model.

    Returns
    -------
    dict of {str: str}
        The compartment whose metabolites each unit may exchange.

    Notes
    -----
    Worked out per unit rather than once for the community, because the
    units need not agree: one model may call it ``e`` and another
    ``extracellular``. cobrapy's own heuristic does the work -- a
    compartment with a recognised name, else the one carrying the most
    boundary reactions -- and ``"e"`` is the fallback if it cannot tell,
    since that is what nearly every reconstruction uses.

    """
    from cobra.medium import find_external_compartment

    if isinstance(given, str):
        return {organism: given for organism in organisms}

    resolved: Dict[str, str] = {}
    for model, organism in zip(units, organisms):
        if given is not None and organism in given:
            resolved[organism] = given[organism]
            continue
        try:
            resolved[organism] = find_external_compartment(model)
        except (RuntimeError, KeyError, IndexError) as error:
            resolved[organism] = "e"
            logger.warning(
                "%s: could not tell which compartment is external (%s), so "
                "assuming 'e'. Pass external_compartment to be sure.",
                organism,
                error,
            )
    distinct = sorted(set(resolved.values()))
    if len(distinct) > 1:
        logger.info("external compartments differ between units: %s", distinct)
    return resolved


def _resolve_bounds(
    bounds: Optional[Mapping[str, Tuple[float, float]]],
    key_of_metabolite: Mapping[str, str],
    what: str,
) -> Dict[str, Tuple[float, float]]:
    """Return bounds keyed by pool identity rather than by metabolite.

    Parameters
    ----------
    bounds : mapping, optional
        Bounds keyed by metabolite identifier, as any unit spells it.
    key_of_metabolite : mapping of {str: str}
        The pool identity of each exchangeable metabolite.
    what : str
        The argument's name, for the warning.

    Returns
    -------
    dict of {str: tuple of float}
        The same bounds, keyed so that any unit's spelling resolves.

    Notes
    -----
    Resolved by looking each identifier up among the metabolites that
    actually exist rather than by stripping a suffix off it, so a
    metabolite named in a medium but present in no unit is reported
    instead of silently doing nothing.

    """
    resolved: Dict[str, Tuple[float, float]] = {}
    unmatched = []
    for met_id, (lower, upper) in (bounds or {}).items():
        key = key_of_metabolite.get(met_id)
        if key is None:
            unmatched.append(met_id)
            continue
        resolved[key] = (float(lower), float(upper))
    if unmatched:
        logger.warning(
            "%d of %d metabolites named in %s are in no unit and were "
            "ignored, for example %s",
            len(unmatched),
            len(bounds or {}),
            what,
            sorted(unmatched)[:5],
        )
    return resolved


def _report_sharing(
    pool_members: Mapping[str, Sequence[str]],
    touched: Mapping[Tuple[str, str], set],
    examples: Mapping[str, Sequence[str]],
) -> None:
    """Warn when the units of a pool have nothing in common.

    Parameters
    ----------
    pool_members : mapping
        The units in each pool.
    touched : mapping of {(str, str): set of str}
        Which units reach each pool metabolite.
    examples : mapping of {str: sequence of str}
        A few exchangeable metabolite identifiers per unit.

    Notes
    -----
    The failure this catches is silent and total. Two models can both
    report an ``e`` compartment and still share nothing, because one
    writes ``glc__D_e`` and the other ``glc_D[e]``: every unit then trades
    through a pool of its own and the community never cross-feeds, with
    no error anywhere. Compartment detection cannot fix that -- the
    identifiers belong to different namespaces -- so it is reported, and
    `metabolite_key` is how to resolve it.

    """
    for pool, members in pool_members.items():
        if len(members) < 2:
            continue
        shared = sum(
            1
            for (where, _), units in touched.items()
            if where == pool and len(units) > 1
        )
        if shared:
            logger.info(
                "pool %r: %d metabolites reached by more than one unit", pool, shared
            )
            continue
        logger.warning(
            "pool %r is shared by %d units but not one metabolite in it is "
            "reached by more than one of them, so nothing can cross-feed. "
            "The usual cause is identifiers from different namespaces; "
            "pass metabolite_key to say how they correspond. Examples: %s",
            pool,
            len(members),
            {unit: list(examples.get(unit, ()))[:3] for unit in list(members)[:3]},
        )


def _measure_anchors(
    model: "Model", anchors: Mapping[str, str], recheck_below: float
) -> Dict[str, float]:
    """Return how much flux each unit's anchor can carry.

    Parameters
    ----------
    model : cobra.Model
        The community model, before coupling is imposed.
    anchors : mapping of {str: str}
        Each unit's anchor reaction.
    recheck_below : float
        Anchors reading less than this after the joint solve are measured
        again on their own.

    Returns
    -------
    dict of {str: float}
        The most flux each unit's anchor can carry, or empty if the model
        could not be solved.

    Notes
    -----
    Measured on the community rather than on each unit alone: a unit that
    cannot grow in isolation may grow perfectly well on what the others
    secrete, and reporting it as dead would be a false alarm. Measured
    before coupling, because coupling a unit to an anchor that cannot
    carry flux pins the unit at zero, which would make the reading agree
    with itself.

    One solve maximising every anchor at once lands on an arbitrary point
    of the optimal face, so a unit can read zero merely for having lost a
    tie to another that uses the same nutrient. That reading is a lower
    bound, which is enough for any unit already clear of
    `recheck_below`; the rest are maximised individually, which costs one
    solve per unit but only for the units in question.

    """
    import numpy as np

    with model:
        model.objective = {
            model.reactions.get_by_id(rxn_id): 1.0 for rxn_id in anchors.values()
        }
        value = model.slim_optimize()
        if value is None or np.isnan(value):
            return {}
        capacity = {
            organism: float(model.reactions.get_by_id(rxn_id).flux)
            for organism, rxn_id in anchors.items()
        }
    for organism, rxn_id in anchors.items():
        # Zero is always worth a second look, whatever the threshold: it
        # is the one reading that is reported, and it is also what a unit
        # shows when it merely lost a tie to another on the same nutrient.
        if capacity[organism] > 0.0 and capacity[organism] >= recheck_below:
            continue
        with model:
            model.objective = model.reactions.get_by_id(rxn_id)
            alone = model.slim_optimize()
        if alone is not None and not np.isnan(alone):
            capacity[organism] = max(capacity[organism], float(alone))
    return capacity


def build_community_model(
    models: Sequence["Model"],
    organisms: Optional[Sequence[str]] = None,
    mode: str = POOLED,
    pools: Optional[Mapping[str, Iterable[str]]] = None,
    environment: Optional[Iterable[str]] = None,
    databases: Optional[Mapping[str, "Model"]] = None,
    biomass_reactions: Optional[Dict[str, str]] = None,
    shared_metabolites: Optional[Iterable[str]] = None,
    link_bounds: Optional[Mapping[str, Tuple[float, float]]] = None,
    pool_medium: Optional[Mapping[str, Tuple[float, float]]] = None,
    shared_compartment: str = SHARED_SUFFIX,
    external_compartment: Optional[Union[str, Mapping[str, str]]] = None,
    metabolite_key: Optional[Callable[["Metabolite"], str]] = None,
    couple: bool = True,
    coupling_factor: float = COUPLING_FACTOR,
    coupling_threshold: float = COUPLING_THRESHOLD,
    check: bool = True,
    model_id: str = "community",
) -> CommunityModel:
    """Join models into one network whose units trade through a pool.

    Parameters
    ----------
    models : sequence of cobra.Model
        The units. They are not modified.
    organisms : sequence of str, optional
        A short identifier for each, in the same order (default: each
        model's own id). These become the suffix on every reaction and
        metabolite, so they must be unique and free of ``__``.
    mode : {"pooled", "shared"}, optional
        How to build the pool (default "pooled"). See the module
        docstring; ``"shared"`` is much smaller but allows only one pool
        and no per-unit bounds on what crosses into it.
    pools : mapping of {str: iterable of str}, optional
        The units meeting in each pool (default: one pool holding all of
        them). A pool is named by a string that is not a unit, and a unit
        may belong to more than one.
    environment : iterable of str, optional
        The pools carrying exchange reactions. Required when there is more
        than one pool: a pool that merely links two units must stay
        closed, or they could bypass the rest of the system.
    databases : mapping of {str: cobra.Model}, optional
        A model of everything a unit may draw on, keyed by unit. The unit
        then holds the whole database, bounded by its own model where the
        two agree, and the extra reactions are recorded in
        `database_reactions` as candidates to gap-fill in.
    biomass_reactions : dict of {str: str}, optional
        The anchor reaction of each unit, keyed by unit (default: the
        objective, or a reaction named ``biomass*``). Only needed when
        `couple` is set.
    shared_metabolites : iterable of str, optional
        Metabolites to leave untagged, so that every unit uses the same
        one. Rarely wanted, but a database with a single ``biomass``
        metabolite drained by one community-wide exchange needs it.
    link_bounds : mapping of {str: tuple of float}, optional
        Bounds on what may cross between a unit and a pool, keyed by
        metabolite as the unit models spell it. This is where measured
        uptake and secretion rates go. Only meaningful when `mode` is
        ``"pooled"``, since otherwise there is nothing to bound.
    pool_medium : mapping of {str: tuple of float}, optional
        Bounds on the exchange reactions between a pool and the
        environment, keyed by metabolite as the unit models spell it.
        Anything unnamed keeps its default bounds.
    shared_compartment : str, optional
        The name of the default pool (default ``"u"``).
    external_compartment : str or mapping, optional
        The compartment whose metabolites are exchangeable: one name for
        every unit, or one per unit. Left out, each unit's is worked out
        from the model, which matters because the units need not agree --
        one may call it ``e`` and another ``extracellular``.
    metabolite_key : callable, optional
        Given a metabolite, return the name of the thing it *is*, so that
        two units naming the same compound differently still meet in the
        pool. The default strips the compartment, which resolves
        ``glc_D_e`` and ``glc_D[e]`` to the same ``glc_D`` but cannot
        reconcile identifiers from different namespaces. Supply this when
        joining models from collections that do not share one.
    couple : bool, optional
        Whether to tie each unit's reactions to its own anchor (default
        True). Without this a unit that is not growing can still run its
        metabolism and feed the others.
    coupling_factor : float, optional
        Flux allowed per unit of anchor flux (default 1000, the published
        value for a biomass anchor). An anchor that runs at a different
        scale needs a different factor; `check` reports what the chosen
        one amounts to.
    coupling_threshold : float, optional
        Flux allowed at zero anchor flux (default 0.01).
    check : bool, optional
        Whether to measure each anchor and report a coupling that would
        throttle its unit (default True). Costs one LP.
    model_id : str, optional
        The identifier of the joined model (default "community").

    Returns
    -------
    CommunityModel
        The joined model and the bookkeeping needed to read it.

    Raises
    ------
    SpectraError
        If the unit names are not usable, the pools do not make sense, an
        anchor cannot be identified, or the result is closed.

    """
    from cobra import Metabolite, Model, Reaction

    if mode not in POOL_MODES:
        raise SpectraError(f"mode must be one of {POOL_MODES}, not {mode!r}.")

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

    pool_members = _resolve_pools(organisms, pools, mode, shared_compartment)
    open_pools = _resolve_environment(environment, pool_members)
    pools_of = {
        organism: tuple(p for p, members in pool_members.items() if organism in members)
        for organism in organisms
    }
    if mode == SHARED and link_bounds:
        raise SpectraError(
            f"link_bounds bounds what crosses between a unit and a pool, but "
            f"mode={SHARED!r} has no such reactions to bound. Use "
            f"mode={POOLED!r}, or put the bounds in pool_medium."
        )

    # Fold each unit's database in before anything is tagged, so that the
    # exchange detection and the tagging see one model per unit.
    units: List["Model"] = []
    draft_of: Dict[str, List[str]] = {}
    database_of: Dict[str, List[str]] = {}
    for model, organism in zip(models, organisms):
        database = (databases or {}).get(organism)
        if database is None:
            units.append(model)
            continue
        merged, draft_ids, database_ids = _merge_with_database(
            model, database, organism
        )
        units.append(merged)
        draft_of[organism] = draft_ids
        database_of[organism] = database_ids

    resolved = dict(biomass_reactions or {})
    anchor_sources: Dict[str, str] = {
        organism: "given" for organism in resolved if organism in set(organisms)
    }
    for model, organism in zip(units, organisms):
        if organism not in resolved:
            resolved[organism], anchor_sources[organism] = _find_anchor(model, organism)

    community = Model(model_id)
    community.solver = models[0].problem
    untagged = set(shared_metabolites or ())
    pooled: Dict[Tuple[str, str], "Metabolite"] = {}
    pool_mets: Dict[str, List[str]] = {pool: [] for pool in pool_members}
    reactions_of: Dict[str, Tuple[str, ...]] = {}
    anchor_ids: Dict[str, str] = {}
    transports: List[Tuple["Reaction", str, "Metabolite", "Metabolite"]] = []
    touched: Dict[Tuple[str, str], set] = {}

    external_of = _external_compartments(units, organisms, external_compartment)

    # Worked out before anything is built, so that a medium can be given
    # in whatever identifiers the caller has and still be matched against
    # the metabolites that exist.
    key_of_metabolite: Dict[str, str] = {}
    examples: Dict[str, List[str]] = {}
    for model, organism in zip(units, organisms):
        outside = external_of[organism]
        for met in model.metabolites:
            if met.compartment != outside:
                continue
            key = (
                metabolite_key(met)
                if metabolite_key is not None
                else _strip_compartment(met.id, outside)
            )
            key_of_metabolite[met.id] = key
            examples.setdefault(organism, []).append(met.id)
    link_by_key = _resolve_bounds(link_bounds, key_of_metabolite, "link_bounds")
    medium_by_key = _resolve_bounds(pool_medium, key_of_metabolite, "pool_medium")

    def tag(identifier: str, organism: str) -> str:
        return f"{identifier}{ORGANISM_SEPARATOR}{organism}"

    for model, organism in zip(units, organisms):
        # Only the exchanges facing the environment are replaced. A demand
        # or sink inside the cell is also a boundary reaction by cobrapy's
        # reckoning, and an anchor that only consumes looks exactly like
        # one, so dropping every boundary reaction would delete it.
        outside = external_of[organism]
        exchanges = {
            rxn.id
            for rxn in model.boundary
            if any(met.compartment == outside for met in rxn.metabolites)
        }

        metabolites: Dict[str, "Metabolite"] = {}
        for met in model.metabolites:
            is_pool_side = mode == SHARED and met.compartment == outside
            if is_pool_side:
                # One row for the whole community, found by what the
                # metabolite is rather than by how this unit spells it,
                # so two spellings of glucose meet instead of passing.
                pool = pools_of[organism][0]
                key = key_of_metabolite[met.id]
                existing = pooled.get((pool, key))
                if existing is None:
                    existing = Metabolite(
                        met.id,
                        formula=met.formula,
                        name=met.name,
                        charge=met.charge,
                        compartment=met.compartment,
                    )
                    community.add_metabolites([existing])
                    pooled[(pool, key)] = existing
                    pool_mets[pool].append(existing.id)
                touched.setdefault((pool, key), set()).add(organism)
                metabolites[met.id] = existing
                continue
            if met.id in untagged:
                existing = (
                    community.metabolites.get_by_id(met.id)
                    if met.id in community.metabolites
                    else None
                )
                if existing is None:
                    existing = Metabolite(
                        met.id,
                        formula=met.formula,
                        name=met.name,
                        charge=met.charge,
                        compartment=met.compartment,
                    )
                    community.add_metabolites([existing])
                metabolites[met.id] = existing
                continue
            copied = Metabolite(
                tag(met.id, organism),
                formula=met.formula,
                name=met.name,
                charge=met.charge,
                compartment=f"{met.compartment}{ORGANISM_SEPARATOR}{organism}",
            )
            metabolites[met.id] = copied
            community.add_metabolites([copied])

        own: List[Tuple["Reaction", "Reaction"]] = []
        for reaction in model.reactions:
            if reaction.id in exchanges:
                # Replaced by a transport into the pool, or by one shared
                # exchange on the pool itself.
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

        if mode == POOLED:
            # One transport per exchangeable metabolite, per pool the unit
            # belongs to: unit's own external compartment <-> pool.
            for met in model.metabolites:
                if met.compartment != outside:
                    continue
                base = key_of_metabolite[met.id]
                for pool in pools_of[organism]:
                    touched.setdefault((pool, base), set()).add(organism)
                    if (pool, base) not in pooled:
                        shared_met = Metabolite(
                            f"{base}_{pool}",
                            formula=met.formula,
                            name=met.name,
                            compartment=pool,
                            charge=met.charge,
                        )
                        pooled[(pool, base)] = shared_met
                        pool_mets[pool].append(shared_met.id)
                        community.add_metabolites([shared_met])
                    shared_met = pooled[(pool, base)]
                    lower, upper = link_by_key.get(base, (-1000.0, 1000.0))
                    transport = Reaction(
                        tag(f"TR_{shared_met.id}", organism),
                        name=f"{met.name} transport, {organism} to {pool}",
                        lower_bound=lower,
                        upper_bound=upper,
                    )
                    transports.append(
                        (transport, organism, metabolites[met.id], shared_met)
                    )

        reactions_of[organism] = tuple(c.id for c, _ in own)
        anchor_ids[organism] = tag(resolved[organism], organism)

    if transports:
        community.add_reactions([t for t, _, _, _ in transports])
        for transport, organism, local, shared_met in transports:
            transport.add_metabolites({local: -1.0, shared_met: 1.0})
            reactions_of[organism] = reactions_of[organism] + (transport.id,)

    # The environment, reachable only through an open pool.
    key_of_pool_metabolite = {
        (pool, met.id): key for (pool, key), met in pooled.items()
    }
    community_exchanges: List[str] = []
    for pool in open_pools:
        for met_id in pool_mets[pool]:
            shared_met = community.metabolites.get_by_id(met_id)
            base = key_of_pool_metabolite[(pool, met_id)]
            lower, upper = medium_by_key.get(base, (0.0, 1000.0))
            exchange = Reaction(
                f"EX_{shared_met.id}",
                name=f"{shared_met.name} community exchange",
                lower_bound=lower,
                upper_bound=upper,
            )
            community.add_reactions([exchange])
            exchange.add_metabolites({shared_met: -1.0})
            community_exchanges.append(exchange.id)

    _report_sharing(pool_members, touched, examples)

    missing = [o for o in organisms if anchor_ids[o] not in community.reactions]
    if missing:
        raise SpectraError(
            f"The anchor reaction of {missing[:5]} is not in the community; "
            f"it may have been an exchange reaction, which is not kept."
        )
    if not community.boundary:
        raise SpectraError(
            "The community has no boundary reaction at all, so nothing can "
            "enter or leave it and no unit can grow. Name a pool in "
            "environment, or keep a demand or sink reaction inside a unit."
        )

    community.objective = {
        community.reactions.get_by_id(anchor_ids[o]): 1.0 for o in organisms
    }

    capacity: Dict[str, float] = {}
    # Nothing to measure until something can be taken up. A community
    # built before its medium is applied is starved by construction, and
    # reporting every unit as inert would be noise rather than news.
    fed = any(rxn.lower_bound < 0 for rxn in community.boundary)
    if check and not fed:
        logger.debug(
            "no boundary reaction permits uptake, so the anchors were not "
            "measured; apply a medium and call coupling_summary() to see them"
        )
    if check and fed:
        # The anchor flux at which coupling stops biting. Anything at or
        # above it needs no second look.
        clear_at = (
            max(0.0, (1.0 - coupling_threshold) / coupling_factor)
            if couple and coupling_factor > 0
            else 0.0
        )
        capacity = _measure_anchors(community, anchor_ids, clear_at)
        for organism in organisms:
            value = capacity.get(organism)
            if value is None:
                continue
            if abs(value) <= 1e-9:
                logger.warning(
                    "%s: its anchor %s cannot carry any flux even when "
                    "maximised on its own, so the unit is inert%s",
                    organism,
                    anchor_ids[organism],
                    " and, coupled, is pinned at zero" if couple else "",
                )
            elif couple and coupling_factor * value + coupling_threshold < 1.0:
                logger.warning(
                    "%s: anchor %s reaches at most %.3g, so coupling caps "
                    "every reaction in %s at %.3g -- raise coupling_factor "
                    "or set couple=False",
                    organism,
                    anchor_ids[organism],
                    value,
                    organism,
                    coupling_factor * value + coupling_threshold,
                )

    if couple:
        constraints = []
        for organism in organisms:
            logger.info(
                "coupling %s on %s (%s); c=%g, u=%g",
                organism,
                anchor_ids[organism],
                anchor_sources.get(organism, "given"),
                coupling_factor,
                coupling_threshold,
            )
            constraints.extend(
                _couple_to_anchor(
                    community,
                    reactions_of[organism],
                    anchor_ids[organism],
                    coupling_factor,
                    coupling_threshold,
                )
            )
        community.add_cons_vars(constraints)
        # Reversed so the context rollback deletes the tail of optlang's variable
        # container first: removing from the front re-indexes everything after it,
        # which is quadratic in the number of auxiliary variables.
        constraints.reverse()

    logger.info(
        "community of %d units in %s mode: %d reactions, %d metabolites, "
        "%d pools, %d community exchanges",
        len(organisms),
        mode,
        len(community.reactions),
        len(community.metabolites),
        len(pool_members),
        len(community_exchanges),
    )
    return CommunityModel(
        model=community,
        organisms=tuple(organisms),
        biomass_reactions=anchor_ids,
        reactions_of=reactions_of,
        community_exchanges=tuple(community_exchanges),
        mode=mode,
        pools={pool: tuple(members) for pool, members in pool_members.items()},
        pool_metabolites={pool: tuple(mets) for pool, mets in pool_mets.items()},
        anchor_sources=anchor_sources,
        anchor_capacity=capacity,
        # Filtered to what survived the build: a database's exchange
        # reactions are replaced by transports, so listing them as
        # candidates would name reactions the model does not have.
        draft_reactions={
            o: tuple(r for r in (tag(i, o) for i in ids) if r in set(reactions_of[o]))
            for o, ids in draft_of.items()
        },
        database_reactions={
            o: tuple(r for r in (tag(i, o) for i in ids) if r in set(reactions_of[o]))
            for o, ids in database_of.items()
        },
        coupling_factor=coupling_factor if couple else None,
        coupling_threshold=coupling_threshold if couple else None,
    )


def build_multi_tissue_model(
    model: "Model",
    tissues: Sequence[str],
    pools: Mapping[str, Iterable[str]],
    environment: Iterable[str],
    anchor_reactions: Optional[Mapping[str, str]] = None,
    **kwargs,
) -> CommunityModel:
    """Replicate one model across tissues and join them.

    Parameters
    ----------
    model : cobra.Model
        The reference network every tissue starts from, typically a
        human genome-scale reconstruction.
    tissues : sequence of str
        The tissue names, which become the tag on every reaction.
    pools : mapping of {str: iterable of str}
        The tissues meeting in each pool. A pool with no tissue of its own
        -- blood, say -- is simply a name that appears here and not in
        `tissues`.
    environment : iterable of str
        The pools that reach the environment. Usually just the blood.
    anchor_reactions : mapping of {str: str}, optional
        The reaction each tissue's coupling gates on, keyed by tissue, as
        named in `model`. A tissue is rarely growing, so this is usually a
        maintenance or demand reaction rather than biomass, and it has to
        be said rather than guessed.
    **kwargs
        Passed to :func:`build_community_model`.

    Returns
    -------
    CommunityModel
        The joined multi-tissue model.

    Notes
    -----
    Every tissue is a copy of the same network, so what distinguishes them
    is their bounds, their core reactions and which pool they can reach.
    The copies are made here rather than by the caller because the
    reference model is usually large and copying it correctly matters.

    """
    # The same model, handed over once per tissue: the builder reads its
    # inputs and never modifies them, so copying a human reconstruction
    # five times would only cost time.
    copies = [model] * len(tissues)
    anchors = (
        {tissue: anchor_reactions[tissue] for tissue in tissues}
        if anchor_reactions
        else None
    )
    return build_community_model(
        copies,
        organisms=list(tissues),
        pools=pools,
        environment=environment,
        biomass_reactions=anchors,
        **kwargs,
    )
