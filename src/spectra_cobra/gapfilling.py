"""Filling the gaps in a draft reconstruction.

A draft usually cannot do everything the organism can: pathways are missing,
so the model fails to grow on a substrate the organism grows on, or cannot
perform a reaction sequence it certainly performs. Gap-filling borrows the
smallest useful set of reactions from a universal model to close that gap.

Two requirements are supported:

* :func:`gapfill_for_growth` requires the model to grow in each of a panel of
  media via a lower bound on the biomass reaction.
* :func:`gapfill_for_tasks` requires the model to perform a list of
  :mod:`~spectra_cobra.tasks`, using :func:`~spectra_cobra.spectra_ccme` in
  Step 1 across any of the four ``problem_type`` formulations
  (``"minNetMILP"``, ``"minNetLP"``, ``"tradeOff"``, or ``"optimBiomass"``)
  followed by a task-by-task ``"minNetMILP"`` repair pass in Step 2.
"""

import math
import warnings
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

import numpy as np
from cobra.util.solver import linear_reaction_coefficients

from ._copy import subset_model
from ._orientation import STOICHIOMETRY
from .exceptions import SpectraError, SpectraSolverError, SpectraValueError
from .extraction import (
    MIN_NET_LP,
    MIN_NET_MILP,
    OPTIM_BIOMASS,
    PROBLEM_TYPES,
    TRADE_OFF,
    _validate_problem_type,
    spectra_ccme,
    spectra_me,
)
from .formulations import (
    DEFAULT_SOLVER_TOLERANCE,
    _validate_positive_finite,
    min_net_milp,
)

if TYPE_CHECKING:
    from cobra.core import Model, Reaction

    from .community import CommunityModel
    from .tasks import MetabolicTask


logger = getLogger(__name__)

#: The growth rate a gap-filled model is asked for when the caller does not
#: say.
DEFAULT_MIN_GROWTH = 0.1

#: Formulations that can be used to gap-fill towards growth.
GROWTH_PROBLEM_TYPES = (MIN_NET_MILP, MIN_NET_LP)

#: What a medium says about one exchange reaction.
Bound = Union[float, Tuple[float, float]]

#: A medium, keyed by exchange reaction.
Medium = Mapping[str, Bound]


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
    skipped : tuple of str
        Tasks the universal model itself cannot perform, and so were skipped.

    """

    model: "Model"
    added: Tuple[str, ...]
    satisfied: Tuple[str, ...] = ()
    unsatisfied: Tuple[str, ...] = ()
    failed: Dict[str, str] = None  # type: ignore[assignment]
    skipped: Tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Default the failure map without making the class mutable."""
        if self.failed is None:
            object.__setattr__(self, "failed", {})

    def summary(self) -> str:
        """Return a one-line account of the result."""
        return (
            f"{len(self.model.reactions)} reactions "
            f"({len(self.added)} added); "
            f"{len(self.satisfied)} satisfied, "
            f"{len(self.unsatisfied)} not, "
            f"{len(self.failed)} unsolved"
        )


def _as_labelled(
    media: Union[Medium, Sequence[Medium], Mapping[str, Medium]],
) -> List[Tuple[str, Medium]]:
    """Return the media as (label, medium) pairs."""
    if isinstance(media, Mapping):
        if media and all(isinstance(v, Mapping) for v in media.values()):
            return list(media.items())
        return [("medium_0", media)]
    return [(f"medium_{index}", medium) for index, medium in enumerate(media)]


def _as_bounds(
    rxn_id: str, bound: "Bound"
) -> Tuple[float, Optional[float]]:
    """Return one medium entry as a lower bound and an optional upper."""
    if isinstance(bound, (int, float)) and not isinstance(bound, bool):
        return float(bound), None
    try:
        lower, upper = bound
        lower, upper = float(lower), float(upper)
    except (TypeError, ValueError) as error:
        raise SpectraValueError(
            f"The medium entry for {rxn_id!r} is {bound!r}, which is neither "
            f"a lower bound nor a (lower, upper) pair."
        ) from error
    if lower > upper:
        raise SpectraValueError(
            f"The medium gives {rxn_id!r} bounds ({lower:g}, {upper:g}), "
            f"which are the wrong way round. Uptake is negative flux, so an "
            f"uptake of at most 10 is a lower bound of -10, and requiring "
            f"between 5 and 10 is (-10, -5)."
        )
    return lower, upper


def _apply_medium(model: "Model", medium: Medium) -> None:
    """Close every exchange, then open the ones the medium names."""
    unknown = [r for r in medium if r not in model.reactions]
    if unknown:
        raise SpectraValueError(
            f"The medium names {len(unknown)} reaction(s) the model does not "
            f"have, for example {sorted(unknown)[:5]}."
        )
    for reaction in model.boundary:
        reaction.lower_bound = 0.0
    for rxn_id, bound in medium.items():
        lower, upper = _as_bounds(rxn_id, bound)
        reaction = model.reactions.get_by_id(rxn_id)
        if upper is None:
            reaction.lower_bound = lower
        else:
            reaction.bounds = (lower, upper)


def _objective_reaction(model: "Model", given: Optional[str]) -> str:
    """Return the reaction whose flux counts as growth."""
    if given is not None:
        if given not in model.reactions:
            raise SpectraValueError(f"No reaction {given!r} in the model.")
        return given
    objective = [r.id for r in model.reactions if r.objective_coefficient != 0]
    if len(objective) == 1:
        return objective[0]
    raise SpectraError(
        f"The model has {len(objective)} reactions in its objective, so which "
        f"one is growth is ambiguous. Pass objective_reaction."
    )


def _grows(model: "Model", medium: Medium, threshold: float) -> bool:
    """Return whether the model reaches the growth threshold in a medium."""
    with model:
        _apply_medium(model, medium)
        value = model.slim_optimize()
    return value is not None and not np.isnan(value) and value >= threshold


def _rebuild(universal: "Model", keep: Set[str]) -> "Model":
    """Return the universal model cut down to the given reactions."""
    return subset_model(universal, set(keep), keep_orphan_metabolites=True)


def _validate_user_weights(
    weights: Optional[Mapping[str, float]],
    univ_rxn_ids: Set[str],
    problem_type: str,
) -> None:
    """Validate user weight overrides against the universal model and `problem_type`."""
    if weights is None:
        return
    if not isinstance(weights, Mapping):
        raise SpectraValueError("weights must be a mapping from reaction ID to weight.")

    allow_negative = problem_type == TRADE_OFF
    for rxn_id, raw_val in weights.items():
        if not isinstance(rxn_id, str) or rxn_id not in univ_rxn_ids:
            raise SpectraValueError(
                f"No reaction {rxn_id!r} in the universal model, so it cannot be weighted."
            )
        if isinstance(raw_val, bool) or not isinstance(raw_val, (int, float, np.number)):
            raise SpectraValueError(
                f"Weight for reaction {rxn_id!r} must be a finite numeric scalar."
            )
        val = float(raw_val)
        if not math.isfinite(val):
            raise SpectraValueError(
                f"Weight for reaction {rxn_id!r} must be a finite numeric scalar."
            )
        if not allow_negative and val < 0.0:
            raise SpectraValueError(
                f"Reaction weights must be non-negative for problem_type "
                f"{problem_type!r} (got {val:g} on {rxn_id!r}); use "
                f"'tradeOff' for signed weights."
            )


def _inclusion_weights(
    reactions: Iterable["Reaction"],
    free: Set[str],
    given: Optional[Mapping[str, float]],
) -> Dict[str, float]:
    """What it costs `gapfill_for_growth` to include each reaction."""
    given = given or {}
    return {
        reaction.id: (
            0.0 if reaction.id in free else float(given.get(reaction.id, 1.0))
        )
        for reaction in reactions
    }


def _build_task_core_weights(
    universal: "Model",
    free_ids: Set[str],
    weights: Optional[Mapping[str, float]],
    problem_type: str,
) -> Dict[str, float]:
    """Build formulation-aware Step-1 weights for `gapfill_for_tasks`."""
    rxn_ids = [rxn.id for rxn in universal.reactions]
    if problem_type == TRADE_OFF:
        # Candidate additions default to -1 penalty; draft/free default to 0.
        core_weights = {
            rid: (0.0 if rid in free_ids else -1.0) for rid in rxn_ids
        }
    elif problem_type == OPTIM_BIOMASS:
        # Keep positive default reward (1.0) on objective reactions even when in draft/free.
        obj_ids = {
            rxn.id
            for rxn, coef in linear_reaction_coefficients(universal).items()
            if coef
        }
        if not obj_ids:
            msg = (
                "universal.objective has no non-zero objective entries; "
                "'optimBiomass' reduces to L1 flux minimization."
            )
            warnings.warn(msg, UserWarning, stacklevel=3)
            logger.warning(msg)
        core_weights = {
            rid: (0.0 if (rid in free_ids and rid not in obj_ids) else 1.0)
            for rid in rxn_ids
        }
    else:
        core_weights = {
            rid: (0.0 if rid in free_ids else 1.0) for rid in rxn_ids
        }

    for rxn_id, val in (weights or {}).items():
        if rxn_id in core_weights:
            core_weights[rxn_id] = float(val)
    return core_weights


def _build_task_repair_costs(
    constrained: "Model",
    universal: "Model",
    already_there: Set[str],
    is_free_now: Set[str],
    weights: Optional[Mapping[str, float]],
    problem_type: str,
) -> Dict[str, float]:
    """Build non-negative Step-2 `minNetMILP` repair costs for `gapfill_for_tasks`."""
    rxn_ids = [rxn.id for rxn in constrained.reactions]
    given = weights or {}

    if problem_type == TRADE_OFF:
        signed_w = {
            rid: (0.0 if rid in is_free_now else -1.0) for rid in rxn_ids
        }
        for rid, val in given.items():
            if rid in signed_w:
                signed_w[rid] = float(val)

        neg_magnitudes = [-w for w in signed_w.values() if w < 0.0]
        min_neg = min([1.0] + neg_magnitudes) if neg_magnitudes else 1.0
        eps_base = 1e-4 * min_neg

        repair_cost: Dict[str, float] = {}
        for rid in rxn_ids:
            w = signed_w[rid]
            if w < 0.0:
                repair_cost[rid] = -w
            elif rid in already_there:
                repair_cost[rid] = 0.0
            else:
                repair_cost[rid] = eps_base / (1.0 + w)
        return repair_cost

    if problem_type == OPTIM_BIOMASS:
        obj_ids = {
            rxn.id
            for rxn, coef in linear_reaction_coefficients(universal).items()
            if coef
        }
        repair_cost = {
            rid: (0.0 if rid in is_free_now else 1.0) for rid in rxn_ids
        }
        for rid, val in given.items():
            if rid not in obj_ids and rid in repair_cost:
                repair_cost[rid] = float(val)
        return repair_cost

    repair_cost = {
        rid: (0.0 if rid in is_free_now else 1.0) for rid in rxn_ids
    }
    for rid, val in given.items():
        if rid in repair_cost:
            repair_cost[rid] = float(val)
    return repair_cost


def gapfill_for_growth(
    draft: "Model",
    universal: "Model",
    media: Union[Medium, Sequence[Medium], Mapping[str, Medium]],
    min_growth: Optional[float] = DEFAULT_MIN_GROWTH,
    objective_reaction: Optional[str] = None,
    tol: float = 1e-4,
    problem_type: str = MIN_NET_MILP,
    free_reactions: Optional[Iterable[str]] = None,
    weights: Optional[Mapping[str, float]] = None,
    time_limit: Optional[float] = 300.0,
    seed: Optional[int] = None,
) -> GapfillResult:
    """Add reactions until a draft grows in every medium.

    Parameters
    ----------
    draft : cobra.Model
        The reconstruction with gaps. Its reactions cost nothing to keep.
    universal : cobra.Model
        The model to borrow reactions from.
    media : dict, sequence of dict, or dict of {str: dict}
        The media to grow in.
    min_growth : float, optional
        The growth rate to demand, applied as a lower bound on the
        objective reaction (default 0.1). Pass None to require only positive
        growth.
    objective_reaction : str, optional
        The reaction whose flux is growth (default: the model's objective,
        if there is exactly one).
    tol : float, optional
        The minimum flux a core reaction must carry (default 1e-4).
    problem_type : {"minNetMILP", "minNetLP"}, optional
        The formulation (default "minNetMILP").
    free_reactions : iterable of str, optional
        Reactions beyond the draft that cost nothing to include.
    weights : mapping of {str: float}, optional
        Non-negative cost weights for candidate reactions.
    time_limit : float, optional
        Seconds per solve (default 300).
    seed : int, optional
        A seed for the extraction's randomised coefficients.

    Returns
    -------
    GapfillResult
        The gap-filled model, what it added, and which media it grows in.

    """
    tol_val = _validate_positive_finite(tol, "tol")
    if problem_type not in GROWTH_PROBLEM_TYPES:
        raise SpectraValueError(
            f"problem_type must be one of {GROWTH_PROBLEM_TYPES} to gap-fill "
            f"towards growth, not {problem_type!r}. tradeOff requires every "
            f"included reaction to carry at least tol, which growth cannot "
            f"satisfy."
        )

    univ_rxn_ids = {rxn.id for rxn in universal.reactions}
    _validate_user_weights(weights, univ_rxn_ids, problem_type)

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
            costs = _inclusion_weights(universal.reactions, free | added, weights)
            try:
                filled = spectra_me(
                    universal,
                    [biomass],
                    tol=tol_val,
                    weights=costs,
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
    """Return the growth rate that counts as growing."""
    return 1e-6 if min_growth is None else float(min_growth)


def gapfill_for_tasks(
    draft: "Model",
    universal: "Model",
    tasks: Sequence["MetabolicTask"],
    tol: float = 1e-4,
    problem_type: str = MIN_NET_MILP,
    free_reactions: Optional[Iterable[str]] = None,
    weights: Optional[Mapping[str, float]] = None,
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
        cannot perform itself are skipped and recorded in ``result.skipped``.
    tol : float, optional
        The minimum flux a core reaction must carry (default 1e-4).
    problem_type : {"minNetMILP", "minNetLP", "tradeOff", "optimBiomass"}, optional
        The formulation the Step-1 core solve uses (default "minNetMILP").
        Step 2 task-by-task repair always uses ``"minNetMILP"``.
    free_reactions : iterable of str, optional
        Reactions beyond the draft that cost nothing to include.
    weights : mapping of {str: float}, optional
        Reaction weight overrides keyed by identifier. For ``"minNetMILP"``,
        ``"minNetLP"``, and ``"optimBiomass"``, weights must be non-negative
        (draft/free default to 0, except objective reactions in
        ``"optimBiomass"`` which default to 1; candidate additions default to
        1). For ``"tradeOff"``, any finite signed real weights are accepted
        (draft/free default to 0, candidate additions to -1) and converted
        monotonically into non-negative repair penalties in Step 2.
    essential_reactions : iterable of str, optional
        Pre-computed essential task reactions to use as the Step-1 core set.
    core_solve : bool, optional
        Whether to run the Step-1 core extraction via :func:`spectra_ccme`
        (default True).
    repair : bool, optional
        Whether to follow up task by task with ``"minNetMILP"`` on any task
        still failing after Step 1 (default True).
    time_limit : float, optional
        Seconds per MILP solve (default 300).
    seed : int, optional
        A seed for the extraction's randomised coefficients.

    Returns
    -------
    GapfillResult
        The gap-filled model, what it added, and which tasks it performs.

    Raises
    ------
    SpectraValueError
        If `tol`, `problem_type`, or `weights` is invalid.

    """
    from .tasks import check_tasks, essential_reactions_for_tasks, task_constraints

    tol_val = _validate_positive_finite(tol, "tol")
    _validate_problem_type(problem_type)

    univ_rxn_ids = {rxn.id for rxn in universal.reactions}
    _validate_user_weights(weights, univ_rxn_ids, problem_type)

    draft_ids = {rxn.id for rxn in draft.reactions if rxn.id in univ_rxn_ids}
    free_ids = set(draft_ids) | (set(free_reactions or ()) & univ_rxn_ids)

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
        core = [r for r in essential_reactions if r in univ_rxn_ids]
    elif core_solve and usable:
        raw_core, _ = essential_reactions_for_tasks(universal, usable)
        core = [r for r in sorted(raw_core) if r in univ_rxn_ids]
    else:
        core = []

    failed: Dict[str, str] = {}
    added: Set[str] = set()

    # Step 1: Core Solve via spectra_ccme so blocked/scaffold reactions never abort Step 1
    if core and core_solve:
        core_weights = _build_task_core_weights(
            universal, free_ids, weights, problem_type
        )
        try:
            filled, _ = spectra_ccme(
                universal,
                core,
                tol=tol_val,
                consistency_type=STOICHIOMETRY,
                weights=core_weights,
                problem_type=problem_type,
                time_limit=time_limit,
                remove_genes=False,
                seed=seed,
            )
            added |= {rxn.id for rxn in filled.reactions} - draft_ids
        except (SpectraSolverError, SpectraError) as error:
            failed["essential_reactions"] = f"{type(error).__name__}: {error}"
            logger.warning("the core solve failed (%s)", error)

    model = _rebuild(universal, draft_ids | added)

    # Step 2: Task-by-task repair pass via minNetMILP with monotone non-negative costs
    if repair and usable:
        free_extra = set(free_reactions or ())
        for result in check_tasks(model, usable):
            if result.ok:
                continue
            present = draft_ids | added
            try:
                with task_constraints(universal, result.task) as constrained:
                    constrained_ids = {r.id for r in constrained.reactions}
                    task_added_rxns = constrained_ids - univ_rxn_ids
                    already_there = present | task_added_rxns
                    is_free_now = already_there | (free_extra & constrained_ids)
                    repair_cost = _build_task_repair_costs(
                        constrained=constrained,
                        universal=universal,
                        already_there=already_there,
                        is_free_now=is_free_now,
                        weights=weights,
                        problem_type=problem_type,
                    )
                    zero_dirs = {r.id: 0 for r in constrained.reactions}
                    solution = min_net_milp(
                        constrained,
                        zero_dirs,
                        repair_cost,
                        tol=tol_val,
                        steady_state=True,
                        time_limit=time_limit,
                        feasibility_tolerance=DEFAULT_SOLVER_TOLERANCE,
                        optimality_tolerance=DEFAULT_SOLVER_TOLERANCE,
                        integrality_tolerance=DEFAULT_SOLVER_TOLERANCE,
                    )
            except (SpectraSolverError, SpectraError) as error:
                failed[result.task.id] = f"{type(error).__name__}: {error}"
                continue
            added |= (solution.included & univ_rxn_ids) - draft_ids
        model = _rebuild(universal, draft_ids | added)

    outcomes = check_tasks(model, usable)
    return GapfillResult(
        model=model,
        added=tuple(sorted(added)),
        satisfied=tuple(r.task.id for r in outcomes if r.ok),
        unsatisfied=tuple(r.task.id for r in outcomes if not r.ok),
        failed=failed,
        skipped=tuple(skipped),
    )


#: What a community exchange costs under ``minNetMILP``.
MILP_EXCHANGE_WEIGHT = 1.0

#: What one costs under ``minNetLP``.
LP_EXCHANGE_WEIGHT = 0.0


def _exchange_weight(given: Optional[float], problem_type: str) -> float:
    """Return what a community exchange should cost."""
    if given is not None:
        return float(given)
    weight = LP_EXCHANGE_WEIGHT if problem_type == MIN_NET_LP else MILP_EXCHANGE_WEIGHT
    logger.info(
        "community exchanges weighted %g, the default for %s", weight, problem_type
    )
    return weight


@dataclass(frozen=True)
class CommunityGapfillResult:
    """What a community gap-fill added, and the model it produced."""

    community: "CommunityModel"
    added: Dict[str, Tuple[str, ...]]
    core: Tuple[str, ...] = ()
    blocked_core: Tuple[str, ...] = ()

    def models(self) -> Dict[str, "Model"]:
        """Return the gap-filled units as models in their own right."""
        return self.community.decompose()

    def summary(self) -> str:
        """Return a one-line account of the result."""
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
    """Return the cost of keeping each reaction of a community."""
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
    """Gap-fill every unit of a community at once."""
    from .community import ORGANISM_SEPARATOR
    from .consistency import consistent_reaction_ids

    model = community.model
    core: Set[str] = set(core_reactions or ())
    if anchors_are_core:
        core |= set(community.biomass_reactions.values())
    unknown = core - {rxn.id for rxn in model.reactions}
    if unknown:
        raise SpectraValueError(
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
