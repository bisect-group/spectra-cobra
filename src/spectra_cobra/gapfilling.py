"""Filling the gaps in a draft reconstruction.

A draft usually cannot do everything the organism can: pathways are missing,
so the model fails to grow on a substrate the organism grows on, or cannot
perform a reaction sequence it certainly performs. Gap-filling borrows the
smallest useful set of reactions from a universal model to close that gap.

Two requirements are supported, and they differ in how the requirement is
expressed rather than in what they do afterwards.

:func:`gapfill_for_growth` requires the model to grow, in each of a panel of
media. Growth is demanded through a *lower bound* on the biomass reaction,
not through `tol`. The distinction matters: `tol` is the flux every core
reaction must carry, so raising it to a growth rate would also demand that
rate of any other core reaction the caller supplies.

:func:`gapfill_for_tasks` requires the model to perform a list of
:mod:`~spectra_cobra.tasks`, by way of the reactions those tasks cannot do
without.
"""

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
    Union,
)

from .exceptions import SpectraError, SpectraSolverError
from .extraction import MIN_NET_LP, MIN_NET_MILP, spectra_me

if TYPE_CHECKING:
    from cobra.core import Model

    from .community import CommunityModel
    from .tasks import MetabolicTask


logger = getLogger(__name__)

#: The growth rate a gap-filled model is asked for when the caller does not
#: say. A real rate rather than a numerical floor: the point of the exercise
#: is a model that grows, and a model scraping along at 1e-6 does not.
#: Replace it with a measured rate where you have one.
DEFAULT_MIN_GROWTH = 0.1

#: Formulations that can be used to gap-fill towards growth. ``tradeOff``
#: is absent on purpose: it requires every included reaction to carry at
#: least `tol`, and growth needs trace fluxes far below that, so it comes
#: back infeasible at every tolerance.
GROWTH_PROBLEM_TYPES = (MIN_NET_MILP, MIN_NET_LP)

#: A medium, as the lower bound to put on each exchange reaction. Anything
#: not named is closed.
Medium = Mapping[str, float]


@dataclass(frozen=True)
class GapfillResult:
    """What a gap-fill added, and the model it produced.

    Parameters
    ----------
    model : cobra.Model
        The draft with the additions applied.
    added : tuple of str
        The reactions taken from the universal model, sorted.
    satisfied : tuple of str
        The requirements the result meets, by label.
    unsatisfied : tuple of str
        The requirements it still does not meet.
    failed : dict of {str: str}
        Requirements whose solve did not finish, and why.

    """

    model: "Model"
    added: Tuple[str, ...]
    satisfied: Tuple[str, ...] = ()
    unsatisfied: Tuple[str, ...] = ()
    failed: Dict[str, str] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        """Default the failure map without making the class mutable."""
        if self.failed is None:
            object.__setattr__(self, "failed", {})

    def summary(self) -> str:
        """Return a one-line account of the result.

        Returns
        -------
        str
            The sizes and outcomes.

        """
        return (
            f"{len(self.model.reactions)} reactions "
            f"({len(self.added)} added); "
            f"{len(self.satisfied)} satisfied, "
            f"{len(self.unsatisfied)} not, "
            f"{len(self.failed)} unsolved"
        )


def _as_labelled(
    media: Union[Sequence[Medium], Mapping[str, Medium]],
) -> List[Tuple[str, Medium]]:
    """Return the media as (label, medium) pairs.

    Parameters
    ----------
    media : sequence of dict, or dict of {str: dict}
        The media, with or without labels.

    Returns
    -------
    list of (str, dict)
        The labelled media.

    """
    if isinstance(media, Mapping):
        return list(media.items())
    return [(f"medium_{index}", medium) for index, medium in enumerate(media)]


def _apply_medium(model: "Model", medium: Medium) -> None:
    """Close every exchange, then open the ones the medium names.

    Parameters
    ----------
    model : cobra.Model
        The model to constrain, inside a ``with`` block.
    medium : dict of {str: float}
        The lower bound for each exchange reaction to open.

    Raises
    ------
    SpectraError
        If the medium names a reaction the model does not have.

    """
    unknown = [r for r in medium if r not in model.reactions]
    if unknown:
        raise SpectraError(
            f"The medium names {len(unknown)} reaction(s) the model does not "
            f"have, for example {sorted(unknown)[:5]}."
        )
    for reaction in model.boundary:
        reaction.lower_bound = 0.0
    for rxn_id, lower in medium.items():
        model.reactions.get_by_id(rxn_id).lower_bound = float(lower)


def _objective_reaction(model: "Model", given: Optional[str]) -> str:
    """Return the reaction whose flux counts as growth.

    Parameters
    ----------
    model : cobra.Model
        The model to inspect.
    given : str, optional
        The caller's choice, if any.

    Returns
    -------
    str
        The reaction identifier.

    Raises
    ------
    SpectraError
        If it is absent, or cannot be identified.

    """
    if given is not None:
        if given not in model.reactions:
            raise SpectraError(f"No reaction {given!r} in the model.")
        return given
    objective = [r.id for r in model.reactions if r.objective_coefficient != 0]
    if len(objective) == 1:
        return objective[0]
    raise SpectraError(
        f"The model has {len(objective)} reactions in its objective, so which "
        f"one is growth is ambiguous. Pass objective_reaction."
    )


def _grows(model: "Model", medium: Medium, threshold: float) -> bool:
    """Return whether the model reaches the growth threshold in a medium.

    Parameters
    ----------
    model : cobra.Model
        The model to test.
    medium : dict of {str: float}
        The medium to apply.
    threshold : float
        The growth rate to reach.

    Returns
    -------
    bool
        Whether it grows that fast.

    """
    import numpy as np

    with model:
        _apply_medium(model, medium)
        value = model.slim_optimize()
    return value is not None and not np.isnan(value) and value >= threshold


def _rebuild(universal: "Model", keep: Set[str]) -> "Model":
    """Return the universal model cut down to the given reactions.

    Parameters
    ----------
    universal : cobra.Model
        The model to cut down.
    keep : set of str
        The reactions to retain.

    Returns
    -------
    cobra.Model
        The reduced model.

    """
    model = universal.copy()
    model.remove_reactions(
        sorted({r.id for r in universal.reactions} - keep), remove_orphans=False
    )
    return model


def gapfill_for_growth(
    draft: "Model",
    universal: "Model",
    media: Union[Sequence[Medium], Mapping[str, Medium]],
    min_growth: Optional[float] = DEFAULT_MIN_GROWTH,
    objective_reaction: Optional[str] = None,
    tol: float = 1e-4,
    problem_type: str = MIN_NET_MILP,
    free_reactions: Optional[Iterable[str]] = None,
    time_limit: Optional[float] = 300.0,
    seed: Optional[int] = None,
) -> GapfillResult:
    """Add reactions until a draft grows in every medium.

    Parameters
    ----------
    draft : cobra.Model
        The reconstruction with gaps. Its reactions cost nothing to keep.
    universal : cobra.Model
        The model to borrow reactions from. The draft should be a subset of
        it, in the sense that the identifiers line up.
    media : sequence of dict, or dict of {str: dict}
        The media to grow in, each mapping an exchange reaction to the
        lower bound to give it. Everything unnamed is closed. Pass a
        mapping to label them.
    min_growth : float, optional
        The growth rate to demand, applied as a lower bound on the
        objective reaction (default 0.1). Pass None to leave the bound
        alone and require only that growth is possible at all.
    objective_reaction : str, optional
        The reaction whose flux is growth (default: the model's objective,
        if there is exactly one).
    tol : float, optional
        The minimum flux a core reaction must carry (default 1e-4). This is
        the usual extraction threshold and is deliberately *not* how the
        growth rate is set; see the notes.
    problem_type : {"minNetMILP", "minNetLP"}, optional
        The formulation (default "minNetMILP").
    free_reactions : iterable of str, optional
        Reactions beyond the draft that cost nothing to include, typically
        those known to be essential for growth in some medium.
    time_limit : float, optional
        Seconds per solve (default 300).
    seed : int, optional
        A seed for the extraction's randomised coefficients.

    Returns
    -------
    GapfillResult
        The gap-filled model, what it added, and which media it grows in.

    Raises
    ------
    SpectraError
        If `problem_type` cannot be used for growth, or the objective is
        ambiguous.

    Notes
    -----
    The growth requirement is a bound on the biomass reaction, not a large
    `tol`. They are easy to confuse, because making biomass the core
    reaction and setting ``tol`` to a growth rate does appear to work. It
    also silently demands that same rate of every *other* core reaction,
    which breaks as soon as the caller supplies one. Keeping the two apart
    additionally makes the result insensitive to ``tol``: on iJO1366 the
    same 15 reactions come back at every tolerance from 1e-4 to 1e-7.

    ``tradeOff`` is not offered. It requires every included reaction to
    carry at least ``tol``, while a growth solution needs trace fluxes well
    below that, so it is infeasible at every tolerance.

    """
    if problem_type not in GROWTH_PROBLEM_TYPES:
        raise SpectraError(
            f"problem_type must be one of {GROWTH_PROBLEM_TYPES} to gap-fill "
            f"towards growth, not {problem_type!r}. tradeOff requires every "
            f"included reaction to carry at least tol, which growth cannot "
            f"satisfy."
        )

    biomass = _objective_reaction(universal, objective_reaction)
    labelled = _as_labelled(media)
    draft_ids = {rxn.id for rxn in draft.reactions}
    boundary_ids = {rxn.id for rxn in universal.boundary}
    free = set(draft_ids) | boundary_ids | set(free_reactions or ())

    added: Set[str] = set()
    failed: Dict[str, str] = {}
    for label, medium in labelled:
        if _grows(_rebuild(universal, draft_ids | added), medium, _floor(min_growth)):
            logger.debug("%s: already grows, nothing to fill", label)
            continue
        with universal:
            _apply_medium(universal, medium)
            if min_growth is not None:
                universal.reactions.get_by_id(biomass).lower_bound = float(min_growth)
            weights = {
                rxn.id: (0.0 if rxn.id in free or rxn.id in added else 1.0)
                for rxn in universal.reactions
            }
            try:
                filled = spectra_me(
                    universal,
                    [biomass],
                    tol=tol,
                    weights=weights,
                    problem_type=problem_type,
                    time_limit=time_limit,
                    seed=seed,
                )
            except (SpectraSolverError, SpectraError) as error:
                failed[label] = f"{type(error).__name__}: {error}"
                logger.warning("%s: could not be gap-filled (%s)", label, error)
                continue
        added |= {rxn.id for rxn in filled.reactions} - draft_ids

    model = _rebuild(universal, draft_ids | added)
    threshold = _floor(min_growth)
    satisfied = [
        label for label, medium in labelled if _grows(model, medium, threshold)
    ]
    unsatisfied = [label for label, _ in labelled if label not in set(satisfied)]
    return GapfillResult(
        model=model,
        added=tuple(sorted(added)),
        satisfied=tuple(satisfied),
        unsatisfied=tuple(unsatisfied),
        failed=failed,
    )


def _floor(min_growth: Optional[float]) -> float:
    """Return the growth rate that counts as growing.

    Parameters
    ----------
    min_growth : float, optional
        The demanded rate, or None.

    Returns
    -------
    float
        The rate to test against.

    """
    return 1e-6 if min_growth is None else float(min_growth)


def gapfill_for_tasks(
    draft: "Model",
    universal: "Model",
    tasks: Sequence["MetabolicTask"],
    tol: float = 1e-4,
    problem_type: str = MIN_NET_MILP,
    free_reactions: Optional[Iterable[str]] = None,
    essential_reactions: Optional[Iterable[str]] = None,
    core_solve: bool = True,
    repair: bool = True,
    time_limit: Optional[float] = 300.0,
    seed: Optional[int] = None,
) -> GapfillResult:
    """Add reactions until a draft can perform a list of metabolic tasks.

    Parameters
    ----------
    draft : cobra.Model
        The reconstruction with gaps.
    universal : cobra.Model
        The model to borrow reactions from.
    tasks : sequence of MetabolicTask
        What the result has to be able to do. Tasks the universal model
        cannot perform itself are skipped, since nothing could be added to
        make them work.
    tol : float, optional
        The minimum flux a core reaction must carry (default 1e-4).
    problem_type : str, optional
        The formulation (default "minNetMILP").
    free_reactions : iterable of str, optional
        Reactions beyond the draft that cost nothing to include.
    essential_reactions : iterable of str, optional
        The reactions essential to the tasks, if you have already computed
        them. Finding them is the expensive part -- 394 s on Human-GEM --
        and it depends only on the universal model and the task list, so it
        is worth doing once and reusing.
    core_solve : bool, optional
        Whether to run the bulk solve that forces the essential reactions
        in (default True). Set it False when the model handed in has
        already been extracted with those reactions as its core, which is
        the usual pipeline: the extraction is yours and carries your
        evidence weights, and this is only here to repair what the core
        could not express.
    repair : bool, optional
        Whether to follow up task by task on whatever still fails (default
        True). This is usually necessary, see the notes.
    time_limit : float, optional
        Seconds per solve (default 300).
    seed : int, optional
        A seed for the extraction's randomised coefficients.

    Returns
    -------
    GapfillResult
        The gap-filled model, what it added, and which tasks it performs.

    Notes
    -----
    There are two phases, and they can be used separately.

    The first turns the tasks into a core reaction set, by way of the
    reactions each task becomes infeasible without, and extracts around it.
    That is a necessary condition and not a sufficient one: a task with two
    alternative routes has no essential reactions at all, so forcing the
    set cannot guarantee it survives. On Human-GEM 34 of 45 tasks have none,
    and the core constraint says nothing whatever about them.

    The second repairs what is left, gap-filling each still-failing task
    inside its own constraints. On iJO1366 that was the difference between
    one failing task and none.

    If you are extracting a context-specific model yourself -- with the
    essential reactions as the core and your own evidence as the weights --
    then the first phase has already happened, and repeating it here would
    both recompute the essential reactions and overrule your weights. Pass
    ``core_solve=False`` and the already-computed ``essential_reactions``:

    .. code-block:: python

       core, _ = essential_reactions_for_tasks(universal, tasks)
       extracted = spectra_me(universal, sorted(core), weights=evidence)
       result = gapfill_for_tasks(
           extracted, universal, tasks,
           essential_reactions=core, core_solve=False,
       )

    """
    from .formulations import min_net_milp
    from .tasks import check_tasks, essential_reactions_for_tasks, task_constraints

    draft_ids = {rxn.id for rxn in draft.reactions}
    free = set(draft_ids) | set(free_reactions or ())

    usable = [
        r.task for r in check_tasks(universal, tasks) if r.ok and not r.task.should_fail
    ]
    skipped = [t.id for t in tasks if t not in usable and not t.should_fail]
    if skipped:
        logger.warning(
            "%d task(s) the universal model cannot perform are skipped: %s",
            len(skipped),
            skipped[:5],
        )

    if essential_reactions is not None:
        core = set(essential_reactions)
    elif core_solve:
        core, _ = essential_reactions_for_tasks(universal, usable)
    else:
        core = set()
    failed: Dict[str, str] = {}
    added: Set[str] = set()

    if core and core_solve:
        weights = {
            rxn.id: (0.0 if rxn.id in free else 1.0) for rxn in universal.reactions
        }
        try:
            filled = spectra_me(
                universal,
                sorted(core),
                tol=tol,
                weights=weights,
                problem_type=problem_type,
                time_limit=time_limit,
                seed=seed,
            )
            added |= {rxn.id for rxn in filled.reactions} - draft_ids
        except (SpectraSolverError, SpectraError) as error:
            failed["essential_reactions"] = f"{type(error).__name__}: {error}"
            logger.warning("the core solve failed (%s)", error)

    model = _rebuild(universal, draft_ids | added)
    if repair:
        for result in check_tasks(model, usable):
            if result.ok:
                continue
            present = draft_ids | added
            with task_constraints(universal, result.task) as constrained:
                weights = {
                    rxn.id: (0.0 if rxn.id in present else 1.0)
                    for rxn in constrained.reactions
                }
                try:
                    solution = min_net_milp(
                        constrained, {}, weights, tol, True, time_limit, None
                    )
                except (SpectraSolverError, SpectraError) as error:
                    failed[result.task.id] = f"{type(error).__name__}: {error}"
                    continue
            added |= {
                r for r in solution.included if r in {x.id for x in universal.reactions}
            } - draft_ids
        model = _rebuild(universal, draft_ids | added)

    outcomes = check_tasks(model, usable)
    return GapfillResult(
        model=model,
        added=tuple(sorted(added)),
        satisfied=tuple(r.task.id for r in outcomes if r.ok),
        unsatisfied=tuple(r.task.id for r in outcomes if not r.ok),
        failed=failed,
    )


#: What a community exchange costs under ``minNetMILP``. Non-zero on
#: purpose: the objective counts reactions, so charging for an exchange
#: makes the gap-filler prefer a neighbour's secretion over fresh material
#: from the medium, which is the point of filling a community at all.
MILP_EXCHANGE_WEIGHT = 1.0

#: What one costs under ``minNetLP``, which is nothing. That objective
#: sums weighted *flux*, and a shared exchange carries the flux of every
#: unit drawing on it, so charging for it penalises a unit for having
#: company. Measured on two hCom organisms: 61 reactions added at weight
#: 1 against 8 at weight 0, where the mixed-integer answer is also 8.
LP_EXCHANGE_WEIGHT = 0.0


def _exchange_weight(given: Optional[float], problem_type: str) -> float:
    """Return what a community exchange should cost.

    Parameters
    ----------
    given : float, optional
        The caller's choice, or None to take the default for the
        formulation.
    problem_type : str
        The formulation being used.

    Returns
    -------
    float
        The weight.

    """
    if given is not None:
        return float(given)
    weight = LP_EXCHANGE_WEIGHT if problem_type == MIN_NET_LP else MILP_EXCHANGE_WEIGHT
    logger.info(
        "community exchanges weighted %g, the default for %s", weight, problem_type
    )
    return weight


@dataclass(frozen=True)
class CommunityGapfillResult:
    """What a community gap-fill added, and the model it produced.

    Parameters
    ----------
    community : CommunityModel
        The gap-filled community, with its bookkeeping narrowed to the
        reactions that survived.
    added : dict of {str: tuple of str}
        The reactions each unit took from its database, named as the
        unit's own model names them.
    core : tuple of str
        The reactions the extraction was required to keep.
    blocked_core : tuple of str
        Core reactions the consistency check found could not carry flux,
        and which were therefore not required.

    """

    community: "CommunityModel"
    added: Dict[str, Tuple[str, ...]]
    core: Tuple[str, ...] = ()
    blocked_core: Tuple[str, ...] = ()

    def models(self) -> Dict[str, "Model"]:
        """Return the gap-filled units as models in their own right.

        Returns
        -------
        dict of {str: cobra.Model}
            One model per unit, untagged, with an exchange wherever the
            unit met a pool.

        """
        return self.community.decompose()

    def summary(self) -> str:
        """Return a one-line account of the result.

        Returns
        -------
        str
            The size of the result and how much was added.

        """
        total = sum(len(ids) for ids in self.added.values())
        return (
            f"{len(self.community.model.reactions)} reactions across "
            f"{len(self.community.organisms)} units; {total} added"
            + (
                f"; {len(self.blocked_core)} core reactions blocked"
                if self.blocked_core
                else ""
            )
        )


def _community_weights(
    community: "CommunityModel",
    exchange_weight: float,
    overrides: Optional[Mapping[str, float]],
) -> Dict[str, float]:
    """Return the cost of keeping each reaction of a community.

    Parameters
    ----------
    community : CommunityModel
        The community to weight.
    exchange_weight : float
        What a community exchange costs.
    overrides : mapping, optional
        Weights to use instead, for any reaction.

    Returns
    -------
    dict of {str: float}
        A weight for every reaction in the model.

    Notes
    -----
    What a unit already has costs nothing, so the gap-fill is free to keep
    all of it; what only its database has costs one, so the solver adds as
    little as it can. Transports between a unit and a pool cost nothing
    either: in the other pooling mode the units simply share the
    compartment, and charging for cross-feeding in one mode but not the
    other would make the two disagree.

    """
    free: Set[str] = set()
    for ids in community.draft_reactions.values():
        free |= set(ids)
    transports = {
        rxn_id
        for unit in community.organisms
        for rxn_id in community.reactions_of[unit]
        if rxn_id.startswith("TR_")
    }
    exchanges = set(community.community_exchanges)

    weights: Dict[str, float] = {}
    for reaction in community.model.reactions:
        if reaction.id in exchanges:
            weights[reaction.id] = float(exchange_weight)
        elif reaction.id in free or reaction.id in transports:
            weights[reaction.id] = 0.0
        else:
            weights[reaction.id] = 1.0
    weights.update({k: float(v) for k, v in (overrides or {}).items()})
    return weights


def gapfill_community(
    community: "CommunityModel",
    core_reactions: Optional[Iterable[str]] = None,
    anchors_are_core: bool = True,
    weights: Optional[Mapping[str, float]] = None,
    exchange_weight: Optional[float] = None,
    tol: float = 1e-4,
    problem_type: str = MIN_NET_MILP,
    consistency_check: bool = True,
    keep_draft: bool = True,
    time_limit: Optional[float] = 300.0,
    seed: Optional[int] = None,
) -> CommunityGapfillResult:
    """Gap-fill every unit of a community at once.

    Parameters
    ----------
    community : CommunityModel
        A community whose units each hold a database of reactions they
        may draw on, as built by passing `databases` to
        :func:`~spectra_cobra.build_community_model`.
    core_reactions : iterable of str, optional
        Reactions the result must keep, named as the community model
        names them.
    anchors_are_core : bool, optional
        Whether each unit's anchor reaction is required (default True).
        This is what makes every unit grow rather than only the ones that
        happen to be cheapest.
    weights : mapping of {str: float}, optional
        Weights to use instead of the defaults, for any reaction. The
        usual reason is sequence evidence: a reaction the organism's
        genome supports should cost less than one it does not.
    exchange_weight : float, optional
        What a community exchange costs. Left out, it depends on the
        formulation: 1 for ``minNetMILP``, which counts reactions, and 0
        for ``minNetLP``, which sums flux and would otherwise charge a
        unit for the uptake of everything sharing the pool with it.
    tol : float, optional
        The minimum flux a core reaction must carry (default 1e-4).
    problem_type : {"minNetMILP", "minNetLP"}, optional
        The formulation (default "minNetMILP"). The LP minimises total
        weighted flux rather than a count, so it is a relaxation: quick,
        and usually larger.
    consistency_check : bool, optional
        Whether to drop the blocked reactions of the joined model before
        gap-filling (default True). Worth its cost: a database reaction
        that cannot carry flux in this community is not a candidate, and
        removing it shrinks the problem.
    keep_draft : bool, optional
        Whether to keep every reaction a unit already had, whether or not
        the solution uses it (default True). This is what separates
        gap-filling from extraction: the formulation returns the smallest
        network that meets the requirement, which would otherwise throw
        away parts of the draft that carry no flux under this one medium.
    time_limit : float, optional
        Seconds to spend on the solve (default 300).
    seed : int, optional
        A seed for the extraction's randomised coefficients.

    Returns
    -------
    CommunityGapfillResult
        The gap-filled community, and what each unit took.

    Raises
    ------
    SpectraError
        If a named core reaction is not in the community.

    Notes
    -----
    Gap-filling the units together rather than one at a time lets a gap
    in one be closed by another's secretion instead of by a new reaction,
    which is why the community result is the smaller one.

    Derived from the community-scale gap-filling of

        S, P. K., Sridhar, S., Alsmadi, N., Mahadevan, R., and Bhatt,
        N. P. (2026). Generalist method to reconstruct metabolic networks
        from multi-omics data at large-scale. bioRxiv.
        https://doi.org/10.64898/2026.04.02.716249

    whose community formulation in turn follows

        Giannari, D., Ho, C. H., and Mahadevan, R. (2021). A gap-filling
        algorithm for prediction of metabolic interactions in microbial
        communities. *PLOS Computational Biology*, 17(6), e1009060.
        https://doi.org/10.1371/journal.pcbi.1009060

    """
    from .community import ORGANISM_SEPARATOR
    from .consistency import consistent_reaction_ids

    model = community.model
    core: Set[str] = set(core_reactions or ())
    if anchors_are_core:
        core |= set(community.biomass_reactions.values())
    unknown = core - {rxn.id for rxn in model.reactions}
    if unknown:
        raise SpectraError(
            f"These core reactions are not in the community model: "
            f"{sorted(unknown)[:5]}."
        )
    all_weights = _community_weights(
        community, _exchange_weight(exchange_weight, problem_type), weights
    )

    working = community
    blocked_core: Tuple[str, ...] = ()
    if consistency_check:
        keep, _ = consistent_reaction_ids(model, tol=tol)
        keep = set(keep)
        logger.info(
            "consistency check: %d of %d reactions can carry flux",
            len(keep),
            len(model.reactions),
        )
        blocked_core = tuple(sorted(core - keep))
        if blocked_core:
            logger.warning(
                "%d core reactions cannot carry flux in this community and "
                "were dropped, for example %s",
                len(blocked_core),
                blocked_core[:5],
            )
        core &= keep
        working = community.with_model(_rebuild(model, keep))

    candidates = sorted(
        rxn_id
        for rxn_id, weight in all_weights.items()
        if weight != 0.0 and rxn_id in {r.id for r in working.model.reactions}
    )
    logger.info(
        "gap-filling %d units with %s: %d candidate reactions, %d core",
        len(working.organisms),
        problem_type,
        len(candidates),
        len(core),
    )
    filled = spectra_me(
        working.model,
        sorted(core),
        tol=tol,
        weights={
            rxn_id: all_weights[rxn_id]
            for rxn_id in (r.id for r in working.model.reactions)
        },
        problem_type=problem_type,
        time_limit=time_limit,
        seed=seed,
        indicator_reactions=candidates if problem_type == MIN_NET_MILP else None,
    )

    kept = {rxn.id for rxn in filled.reactions}
    if keep_draft:
        present = {rxn.id for rxn in working.model.reactions}
        drafted = {
            rxn_id
            for ids in working.draft_reactions.values()
            for rxn_id in ids
            if rxn_id in present
        }
        logger.info(
            "keeping %d draft reactions the solution does not use",
            len(drafted - kept),
        )
        kept |= drafted
        result = working.with_model(_rebuild(working.model, kept))
    else:
        result = working.with_model(filled)
    added: Dict[str, Tuple[str, ...]] = {}
    for unit in community.organisms:
        suffix = f"{ORGANISM_SEPARATOR}{unit}"
        taken = sorted(
            rxn_id[: -len(suffix)] if rxn_id.endswith(suffix) else rxn_id
            for rxn_id in community.database_reactions.get(unit, ())
            if rxn_id in kept
        )
        added[unit] = tuple(taken)
    return CommunityGapfillResult(
        community=result,
        added=added,
        core=tuple(sorted(core)),
        blocked_core=blocked_core,
    )
