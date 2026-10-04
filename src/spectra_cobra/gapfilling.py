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
