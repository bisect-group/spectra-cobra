"""The smallest network that still does the job.

A genome-scale reconstruction is hugely redundant. Most of its reactions
are not needed to grow on any one carbon source, and a large subnetwork
can be deleted without the model noticing. The *minimal reactome* is the
smallest set of reactions that still meets a stated requirement -- the
original question of Burgard, Vaidyaraman and Maranas (2001), who asked
how few reactions *Escherichia coli* needs to grow, and found the answer
depends entirely on what is in the medium: 224 reactions on glucose
alone, 229 on acetate alone, 122 on a medium engineered to supply as much
as possible.

That is a minimum-cardinality problem, so it is a MILP. Give every
reaction a binary :math:`y_j` tied to its flux by
:math:`lb_j y_j \\le v_j \\le ub_j y_j`, so that :math:`y_j = 0` pins the
flux to zero, and minimise :math:`\\sum_j w_j y_j` subject to
:math:`S v = 0` and whatever the network has to achieve. That is
:func:`~spectra_cobra.min_net_milp` with no core reactions at all: the
requirements enter as bounds rather than as a core set.

What the network has to achieve can be any of

* a **medium** -- which exchanges are open, and how fast,
* a **minimum growth rate** -- a lower bound on the biomass reaction,
* **production** of one or more metabolites, at a stated rate,
* a list of **metabolic tasks** it must still be able to perform.

The first three are one condition and share one flux vector. A task is
not: :func:`~spectra_cobra.tasks.task_constraints` closes the model's own
exchanges and opens its own, so each task is a *different environment*,
and a network that performs two tasks does so with two different flux
distributions. The reactions are shared but the fluxes are not. So with
tasks the problem becomes one flux vector per condition over a single set
of binaries:

.. math::

   \\min_y \\sum_j w_j y_j \\quad\\text{s.t.}\\quad
   S^{(k)} v^{(k)} = 0, \\;
   lb^{(k)}_j y_j \\le v^{(k)}_j \\le ub^{(k)}_j y_j
   \\quad \\forall k

This is the "protected functions" idea of NetworkReducer and MinNW: a
reduced network must keep not one phenotype but a list of them. Solving
the conditions jointly is what makes the answer minimal *across* them.
Filling them one at a time and taking the union, as
:func:`~spectra_cobra.gapfill_for_tasks` does, gives a feasible network
but not the smallest one.

References
----------
Burgard, A. P., Vaidyaraman, S., and Maranas, C. D. (2001). Minimal
reaction sets for *Escherichia coli* metabolism under different growth
requirements and uptake environments. *Biotechnology Progress*, 17(5),
791-797. https://doi.org/10.1021/bp0100880

Sambamoorthy, G., and Raman, K. (2020). MinReact: a systematic approach
for identifying minimal metabolic networks. *Bioinformatics*, 36(15),
4309-4315. https://doi.org/10.1093/bioinformatics/btaa497

Erdrich, P., Steuer, R., and Klamt, S. (2015). An algorithm for the
reduction of genome-scale metabolic network models to meaningful core
models. *BMC Systems Biology*, 9, 48.
https://doi.org/10.1186/s12918-015-0191-x

Röhl, A., and Bockmayr, A. (2017). A mixed-integer linear programming
approach to the reduction of genome-scale metabolic networks. *BMC
Bioinformatics*, 18, 2. https://doi.org/10.1186/s12859-016-1412-z
"""

from contextlib import contextmanager
from dataclasses import dataclass, field
from logging import getLogger
from typing import (
    TYPE_CHECKING,
    Dict,
    Iterable,
    Iterator,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)

from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

from ._copy import subset_model
from .exceptions import SpectraError, SpectraInfeasibleCoreError
from .formulations import INCLUSION_CUTOFF_FACTOR

if TYPE_CHECKING:
    from cobra.core import Model

    from .tasks import MetabolicTask


logger = getLogger(__name__)

#: What separates a reaction's identifier from the condition it belongs to
#: in the replicated model. Only ever seen in solver output, since the
#: result is reported in the model's own identifiers.
CONDITION_SEPARATOR = "__"

#: The label of the condition carrying the medium, the growth requirement
#: and the production requirements. The tasks supply the rest.
GROWTH_CONDITION = "growth"

#: The fraction of what the model can achieve in the medium that the
#: minimal reactome must still achieve, when nothing else is said. One
#: means lose no growth at all, which is MinReact's own default.
DEFAULT_GROWTH_FRACTION = 1.0

#: The absolute flux at which a reaction counts as used, whatever its
#: binary says. Deliberately far below any solver's own tolerance: a
#: binary the solver leaves at its integrality tolerance rather than at
#: exactly zero still multiplies a flux bound of up to a thousand, so a
#: reaction the solution genuinely depends on can be reported as switched
#: off. Dropping it breaks the chain it was part of.
INCLUSION_CUTOFF = 1e-11

#: The flux below which a reaction counts as unused when shrinking the
#: candidates for the essentiality search. Matches the tasks module.
FLUX_EPSILON = 1e-9

#: How far below the requirement a verified result may fall and still
#: count as meeting it, as a fraction. Reading a rate off a second solve
#: of a rebuilt model does not reproduce the first solve bit for bit.
VERIFICATION_SLACK = 1e-6


@dataclass(frozen=True)
class MinimalReactome:
    """The smallest network found, and what it still does.

    Parameters
    ----------
    model : cobra.Model
        The minimal reactome, as a model.
    kept : tuple of str
        The reactions it retained, sorted.
    removed : tuple of str
        The reactions dropped, sorted. Includes the blocked reactions
        preprocessing deleted, which were never candidates.
    growth : dict of {str: float}
        The rate the minimal reactome reaches in each medium, keyed by the
        medium's label. Empty when no growth was required.
    required_growth : dict of {str: float}
        The rate it had to reach in each.
    production : dict of {str: dict of {str: float}}
        The flux each required product carries in each medium, in the same
        solution the growth rate was read from.
    required_production : dict of {str: float}
        The rate each product had to carry, the same in every medium.
    conditions : tuple of str
        The labels of the conditions imposed, in the order solved.
    essential : tuple of str
        The reactions preprocessing proved are in every possible answer,
        because some condition is infeasible without them. Empty when
        `preprocess` was off.
    kept_but_blocked : dict of {str: tuple of str}
        Per growth condition, the `keep_reactions` entries that cannot
        carry `tol` of flux in it. They are still in the result, but
        nothing forces them to do anything there.
    unsatisfied : tuple of str
        The conditions the result does not actually meet. Empty is the
        expected outcome; anything here means the verification caught
        something the solve did not.

    """

    model: "Model"
    kept: Tuple[str, ...]
    removed: Tuple[str, ...]
    growth: Dict[str, float] = field(default_factory=dict)
    required_growth: Dict[str, float] = field(default_factory=dict)
    production: Dict[str, Dict[str, float]] = field(default_factory=dict)
    required_production: Dict[str, float] = field(default_factory=dict)
    conditions: Tuple[str, ...] = ()
    essential: Tuple[str, ...] = ()
    kept_but_blocked: Dict[str, Tuple[str, ...]] = field(default_factory=dict)
    unsatisfied: Tuple[str, ...] = ()

    def summary(self) -> str:
        """Return a one-line account of the result.

        Returns
        -------
        str
            The size, the growth rates and anything unmet.

        """
        total = len(self.kept) + len(self.removed)
        parts = [f"{len(self.kept)} of {total} reactions"]
        for label, reached in sorted(self.growth.items()):
            needed = self.required_growth.get(label, 0.0)
            name = "growth" if len(self.growth) == 1 else f"growth on {label}"
            parts.append(f"{name} {reached:.4g} (needed {needed:.4g})")
        if self.required_production:
            wanted = ", ".join(
                f"{r} >= {v:.4g}" for r, v in sorted(self.required_production.items())
            )
            parts.append(f"producing {wanted}")
        if len(self.conditions) > 1:
            parts.append(f"{len(self.conditions)} conditions")
        if self.essential:
            parts.append(f"{len(self.essential)} essential")
        if self.kept_but_blocked:
            stuck = sorted({r for rs in self.kept_but_blocked.values() for r in rs})
            parts.append(f"{len(stuck)} kept but unable to carry flux")
        if self.unsatisfied:
            parts.append(f"UNMET: {', '.join(self.unsatisfied)}")
        return "; ".join(parts)


Medium = Mapping[str, float]

#: One medium, several unlabelled, or several labelled.
Media = Union[Medium, Sequence[Medium], Mapping[str, Medium]]


def minimal_reactome(
    model: "Model",
    medium: Optional[Media] = None,
    biomass_reaction: Optional[str] = None,
    min_growth: Optional[float] = None,
    growth_fraction: Optional[float] = None,
    products: Optional[Mapping[str, float]] = None,
    tasks: Optional[Sequence["MetabolicTask"]] = None,
    keep_reactions: Optional[Iterable[str]] = None,
    weights: Optional[Mapping[str, float]] = None,
    tol: float = 1e-4,
    preprocess: bool = True,
    n_solutions: int = 1,
    time_limit: Optional[float] = 7200.0,
    inclusion_cutoff: Optional[float] = None,
    seed: Optional[int] = None,
    check: bool = True,
) -> Union[MinimalReactome, List[MinimalReactome]]:
    """Find the smallest set of reactions that still meets a requirement.

    Parameters
    ----------
    model : cobra.Model
        The model to cut down. It is not modified.
    medium : dict, sequence of dict, or dict of {str: dict}, optional
        The exchange reactions to open and the uptake rate to allow each,
        as a lower bound. Everything unnamed is closed. Default None,
        meaning leave the model's own bounds alone.

        Pass **several media** to require that the one network works in all
        of them, as a sequence or as a mapping from label to medium. Each
        becomes its own condition with its own flux vector, exactly as a
        task does, so the result is the smallest network that grows
        aerobically *and* anaerobically rather than the smaller of the two.
        `growth_fraction` is measured separately in each.
    biomass_reaction : str, optional
        The reaction whose flux is growth (default: the model's objective,
        if there is exactly one).
    min_growth : float, optional
        The growth rate to demand, as an absolute lower bound on the
        biomass reaction. Mutually exclusive with `growth_fraction`.
    growth_fraction : float, optional
        The growth rate to demand, as a fraction of what the model
        achieves in each medium before anything is removed (default 1.0,
        i.e. lose no growth at all). Pass 0 to not require growth.
    products : dict of {str: float}, optional
        Metabolites or reactions the minimal reactome must still make, and
        the flux to demand of each, as a lower bound. A reaction identifier
        is bounded directly; a metabolite identifier is resolved to the
        boundary reaction the model already has for it. A metabolite the
        model does not exchange is refused rather than given one -- see the
        notes. These are required *alongside* growth, in the same flux
        distribution.
    tasks : sequence of MetabolicTask, optional
        Functions the minimal reactome must still be able to perform. Each
        becomes its own condition, with its own flux vector and its own
        medium, sharing the one set of binaries. See the notes.
    keep_reactions : iterable of str, optional
        Reactions to retain whatever they do: they are given no binary, so
        they cost nothing and cannot be removed. MinReact's "reactions to
        be retained".
    weights : dict of {str: float}, optional
        What each reaction costs to keep, keyed by identifier; anything
        left out costs 1. Non-negative. A cost of 0 still leaves the
        reaction removable -- use `keep_reactions` to force one in.
    tol : float, optional
        The flux threshold (default 1e-4). Used by the consistency check
        that finds the blocked reactions, and it sets the default
        `inclusion_cutoff`.
    preprocess : bool, optional
        Whether to decide what can be decided before the MILP starts
        (default True). Two passes, both exact: reactions that can carry
        no flux in *any* condition are deleted outright, and reactions
        some condition is infeasible without are proved essential and
        given no binary. On iJO1366 on glucose that is 878 blocked and 407
        essential, which is half the binaries. See the notes.
    seed : int, optional
        A seed for the consistency check's random coefficients, making
        preprocessing reproducible (default None).
    n_solutions : int, optional
        How many alternative minimal reactomes to return (default 1).
        Each is forbidden from recurring exactly before the next solve,
        so they are distinct but not necessarily all of them, and only
        the first is guaranteed minimal.
    time_limit : float, optional
        Seconds to allow each solve (default 7200). A solution found
        before the limit is used even if minimality was not proven.
    inclusion_cutoff : float, optional
        The absolute flux at which a reaction is kept whatever its binary
        says (default ``tol * 1e-7``, so 1e-11 at the default `tol`, which
        is the same rule :func:`~spectra_cobra.spectra_me` uses). See
        :data:`INCLUSION_CUTOFF`; raising it is how you would get a
        strictly minimal count at the risk of a network that cannot carry
        flux.
    check : bool, optional
        Whether to verify the result against every requirement and report
        what it misses (default True).

    Returns
    -------
    MinimalReactome, or list of MinimalReactome
        The minimal reactome, or `n_solutions` of them.

    Raises
    ------
    SpectraError
        If the arguments conflict, name something the model does not have,
        or demand growth the model cannot reach in the medium.

    See Also
    --------
    spectra_cobra.spectra_me : Extraction around a core set of reactions.
    spectra_cobra.minimal_microbiome : The same question over organisms.

    Notes
    -----
    The medium is the whole story. Burgard and colleagues found 224
    reactions sufficient for *E. coli* on glucose and 122 on a medium
    supplying everything it could transport: a minimal reactome is minimal
    *for an environment*, never in the abstract.

    **Nothing is added to the network.** Requiring the production of a
    metabolite the model does not exchange would need a drain for it, and
    a drain added here would make the result the smallest subnetwork of a
    model you do not have -- feasible only because of the reaction that
    was added to measure it. If that is what you want, add the reaction to
    your own model first, where it is visible and counted like any other.

    Preprocessing decides two things exactly, so it cannot change the
    *optimum* -- on iJO1366 on glucose the objective is 28 either way, the
    435-reaction answer being those 28 plus 407 essential. What it can
    change by a reaction or two is the network finally reported, because
    the trace-flux rule behind `inclusion_cutoff` keeps whatever carries
    flux with its binary off, and a differently-shaped relaxation leaves
    traces in different places. A reaction that can carry no
    flux in any condition is in no feasible solution, so it is deleted. A
    reaction some condition is *infeasible without* is in every feasible
    solution, so it is forced in and given no binary -- and because the
    growth condition carries the biomass lower bound, "infeasible without"
    there means "cannot reach the required growth rate without", which is
    the right notion rather than a growth comparison. Both are unions and
    intersections over conditions respectively: essential **anywhere** is
    essential; blocked **everywhere** is blocked.

    Growth, production and the medium are one condition, solved as one
    flux vector. Tasks are not, because a task defines its own medium. The
    conditions are therefore replicated -- one copy of the network each --
    and tied together by one binary per reaction, so what is minimised is
    the size of the union while every condition is satisfied separately.
    The model handed to the solver is `k` times the size of yours, which
    is what makes a long task list expensive.

    This is the opposite operation to :func:`~spectra_cobra.gapfill_for_tasks`
    and they use the same vocabulary of tasks. Gap-filling adds reactions
    until the tasks pass; this removes reactions while they still do.

    Examples
    --------
    How few reactions does the textbook *E. coli* core model need to grow
    at half its maximum on glucose?

    >>> from cobra.io import load_model
    >>> from spectra_cobra import minimal_reactome
    >>> model = load_model("textbook")  # doctest: +SKIP
    >>> result = minimal_reactome(  # doctest: +SKIP
    ...     model,
    ...     medium={"EX_glc__D_e": -10.0, "EX_o2_e": -1000.0,
    ...             "EX_nh4_e": -1000.0, "EX_pi_e": -1000.0,
    ...             "EX_h2o_e": -1000.0, "EX_h_e": -1000.0,
    ...             "EX_co2_e": -1000.0},
    ...     growth_fraction=0.5,
    ... )
    >>> print(result.summary())  # doctest: +SKIP

    """
    if min_growth is not None and growth_fraction is not None:
        raise SpectraError(
            "Pass min_growth or growth_fraction, not both: one is an "
            "absolute rate and the other a fraction of what the model "
            "achieves, so giving both leaves it ambiguous which applies."
        )
    if n_solutions < 1:
        raise SpectraError(f"n_solutions must be at least 1, not {n_solutions}.")

    # A full copy, built rather than cloned: ``subset_model`` keeping
    # everything is faster than ``copy_model`` and gives the same model.
    base = subset_model(
        model,
        {rxn.id for rxn in model.reactions},
        remove_genes=False,
        keep_orphan_metabolites=True,
    )
    original_ids = {rxn.id for rxn in base.reactions}
    tasks = list(tasks or ())
    if inclusion_cutoff is None:
        inclusion_cutoff = tol * INCLUSION_CUTOFF_FACTOR
    # Asking for a growth rate of zero is asking for no growth, so it does
    # not oblige the caller to say which reaction growth would have been.
    growth_requested = (min_growth or 0.0) > 0.0 or (growth_fraction or 0.0) > 0.0
    biomass = _biomass_reaction(base, biomass_reaction, growth_requested)

    labelled_media = _as_labelled_media(medium)
    for _, one in labelled_media:
        if one is not None:
            _check_medium(base, one)

    required_products = _resolve_products(base, products or {})
    required_growth = {
        label: _required_growth(base, one, biomass, min_growth, growth_fraction, label)
        for label, one in labelled_media
    }

    requested = set(keep_reactions or ())
    if requested:
        from .extraction import _warn_on_loose_tolerance

        _warn_on_loose_tolerance(base, tol)
    unknown = requested - original_ids
    if unknown:
        raise SpectraError(
            f"keep_reactions names {len(unknown)} reaction(s) the model does "
            f"not have, for example {sorted(unknown)[:5]}."
        )

    if (
        not any(r > 0.0 for r in required_growth.values())
        and not required_products
        and not tasks
    ):
        raise SpectraError(
            "Nothing was asked of the minimal reactome, so the smallest "
            "network meeting the requirement is the empty one. Give a growth "
            "requirement, a product or a task list. A medium on its own is "
            "not a requirement -- it says what may be taken up, not what has "
            "to happen."
            + (
                " This model has no unambiguous objective reaction, so growth "
                "was not required of it; pass biomass_reaction to say which "
                "reaction growth is."
                if biomass is None
                else ""
            )
        )

    _check_tasks_are_possible(base, tasks)
    for label, one in labelled_media:
        _check_growth_is_possible(
            base, one, biomass, required_growth[label], required_products, label
        )
    forceable, kept_but_blocked = _forceable_keeps(
        base,
        labelled_media,
        biomass,
        required_growth,
        required_products,
        requested,
        tol,
    )
    conditions = _conditions(
        base, labelled_media, biomass, required_growth, required_products, tasks
    )

    blocked, essential = (
        _preprocess(base, conditions, original_ids, requested, tol, seed)
        if preprocess
        else (set(), set())
    )
    if blocked:
        base = subset_model(
            base,
            {rxn.id for rxn in base.reactions} - blocked,
            keep_orphan_metabolites=True,
        )
    protected = (requested | essential) - blocked

    counted = [rxn.id for rxn in base.reactions if rxn.id not in protected]
    costs = _costs(weights, counted, protected)
    joint, groups, copies = _replicate(base, conditions, counted)
    core = sorted(
        copies[label][rxn_id]
        for label, usable in forceable.items()
        if label in copies
        for rxn_id in usable
    )
    logger.info(
        "minimising %d reactions over %d condition(s) (%s): %d reactions and "
        "%d metabolites in the joint problem, %d core",
        len(counted),
        len(conditions),
        ", ".join(label for label, _ in conditions),
        len(joint.reactions),
        len(joint.metabolites),
        len(core),
    )

    from .extraction import MIN_NET_MILP, PATHWAY_EXCLUSION, spectra_me

    try:
        _, solutions = spectra_me(
            joint,
            core,
            tol=tol,
            weights=costs,
            n_solutions=n_solutions,
            alt_solution_method=PATHWAY_EXCLUSION,
            problem_type=MIN_NET_MILP,
            time_limit=time_limit,
            remove_genes=False,
            seed=seed,
            inclusion_cutoff=inclusion_cutoff,
            indicator_groups=groups,
            return_solutions=True,
        )
    except SpectraInfeasibleCoreError as error:
        raise SpectraInfeasibleCoreError(
            f"The reactions in keep_reactions cannot all carry {tol:.3g} of "
            f"flux at once, although each can on its own. {error}"
        ) from error
    if n_solutions == 1:
        solutions = [solutions]

    origin = {
        copy_id: rxn_id
        for per_condition in copies.values()
        for rxn_id, copy_id in per_condition.items()
    }
    results: List[MinimalReactome] = []
    for index, solution in enumerate(solutions):
        kept = {origin[r] for r in solution.included if r in origin} | protected
        results.append(
            _report(
                base,
                original_ids,
                kept,
                labelled_media,
                biomass,
                required_growth,
                required_products,
                tasks,
                tuple(label for label, _ in conditions),
                tuple(sorted(essential)),
                forceable,
                kept_but_blocked,
                tol,
                check,
            )
        )
        logger.info("solution %d: %s", index + 1, results[-1].summary())

    return results[0] if n_solutions == 1 else results


def _as_labelled_media(medium: Optional[Media]) -> List[Tuple[str, Optional[Medium]]]:
    """Return the media as (label, medium) pairs, however they were given.

    Parameters
    ----------
    medium : dict, sequence of dict, or dict of {str: dict}, optional
        One medium, several, or several with labels. None means leave the
        model's own bounds alone.

    Returns
    -------
    list of (str, dict or None)
        One entry per growth condition.

    Raises
    ------
    SpectraError
        If a labelled mapping mixes media with bounds, or a sequence is
        empty.

    Notes
    -----
    A single medium maps exchange identifiers to floats; a labelled set
    maps labels to media. The two are told apart by the type of the
    values, which is unambiguous because a bound is never a mapping.

    """
    if medium is None:
        return [(GROWTH_CONDITION, None)]

    if isinstance(medium, Mapping):
        values = list(medium.values())
        mappings = [isinstance(value, Mapping) for value in values]
        if values and all(mappings):
            return [(str(label), one) for label, one in medium.items()]
        if any(mappings):
            raise SpectraError(
                "The medium mixes uptake bounds with whole media, so whether "
                "it is one medium or several cannot be told. Give either "
                "{exchange: bound} or {label: {exchange: bound}}."
            )
        return [(GROWTH_CONDITION, medium)]

    media = list(medium)
    if not media:
        raise SpectraError(
            "An empty sequence of media was given. Pass None to leave the "
            "model's own bounds alone, or {} for a medium that supplies "
            "nothing."
        )
    if len(media) == 1:
        return [(GROWTH_CONDITION, media[0])]
    return [(f"medium_{index}", one) for index, one in enumerate(media)]


def _preprocess(
    base: "Model",
    conditions: Sequence[Tuple[str, object]],
    original_ids: Set[str],
    requested: Set[str],
    tol: float,
    seed: Optional[int],
) -> Tuple[Set[str], Set[str]]:
    """Decide before the MILP what the MILP would have had to decide.

    Parameters
    ----------
    base : cobra.Model
        The working copy.
    conditions : sequence of (str, callable)
        The conditions, as returned by :func:`_conditions`.
    original_ids : set of str
        The reactions of the model, so that the scaffolding a task adds is
        never reported as blocked or essential.
    requested : set of str
        The caller's `keep_reactions`, which are never worth testing since
        they are kept regardless.
    tol : float
        The flux threshold for the consistency check.
    seed : int, optional
        A seed for its random coefficients.

    Returns
    -------
    set of str
        The reactions blocked in **every** condition, which no feasible
        solution of any condition uses and which can be deleted.
    set of str
        The reactions **some** condition is infeasible without, which every
        feasible solution contains and which need no binary.

    Notes
    -----
    Both passes are exact, so the optimum is unchanged -- only the number
    of binaries the solver has to branch on. Blocked is an intersection
    over conditions and essential a union, and the asymmetry is the point:
    a reaction useless under the growth medium may be the only route a
    task has, while a reaction some task cannot do without has to be in
    the answer whatever the growth condition thinks of it.

    """
    from .consistency import consistent_reaction_ids

    blocked: Optional[Set[str]] = None
    for label, factory in conditions:
        with factory() as constrained:
            # Count a reaction as able to carry flux if it can carry *any*,
            # not if it can reach `tol`. The difference decides whether a
            # reaction carrying a load-bearing trace is deleted here, and
            # deleting one breaks the chain it belonged to -- the same trap
            # the inclusion cutoff exists for. Erring towards keeping costs
            # one binary; erring the other way costs correctness.
            consistent, n_lps = consistent_reaction_ids(
                constrained,
                tol=tol,
                detection_cutoff=constrained.tolerance,
                seed=seed,
            )
        here = original_ids - consistent
        logger.info(
            "%s: %d of %d reactions cannot carry flux [%d LPs]",
            label,
            len(here),
            len(original_ids),
            n_lps,
        )
        blocked = here if blocked is None else (blocked & here)
    blocked = (blocked or set()) - requested

    candidates = original_ids - blocked - requested
    essential: Set[str] = set()
    for label, factory in conditions:
        remaining = candidates - essential
        if not remaining:
            break
        with factory() as constrained:
            found = _essential_in(constrained, remaining)
        logger.info("%s: %d reactions are essential to it", label, len(found))
        essential |= found

    logger.info(
        "preprocessing: %d blocked everywhere (deleted), %d essential "
        "somewhere (no binary), %d binaries left to decide",
        len(blocked),
        len(essential),
        len(original_ids) - len(blocked) - len(essential) - len(requested),
    )
    return blocked, essential


def _essential_in(constrained: "Model", candidates: Set[str]) -> Set[str]:
    """Return the reactions this condition becomes infeasible without.

    Parameters
    ----------
    constrained : cobra.Model
        A model already inside its condition's context, so every bound the
        condition imposes -- the medium, the lower bound on biomass, the
        product floors, or a task's own exchanges -- is in place.
    candidates : set of str
        The reactions to consider.

    Returns
    -------
    set of str
        The identifiers of the essential ones.

    Notes
    -----
    The test is feasibility, not growth: with the requirements already
    expressed as bounds, "the condition is infeasible without this
    reaction" *is* "the required growth rate cannot be reached without
    it". Comparing growth rates would be a second, redundant way of
    saying the same thing, and would not transfer to a task, which has no
    growth at all.

    Only a reaction carrying flux in some feasible solution can be
    essential, so the candidates are first shrunk by minimising total
    flux: anything a parsimonious solution leaves at zero has an
    alternative and is dropped, many at a time, for one LP. Only the
    survivors are knocked out one by one. This is the same strategy
    :func:`~spectra_cobra.essential_reactions_for_task` uses.

    """
    from cobra.flux_analysis.parsimonious import add_pfba

    pool = {r for r in candidates if r in constrained.reactions}
    constrained.objective = constrained.problem.Objective(0.0)
    if not _feasible(constrained):
        return set()

    # A reaction whose own bounds exclude zero cannot be switched off at
    # all, so it is in every solution by construction. Knocking one out
    # would not even be the right question: the biomass reaction carries
    # the growth requirement in its lower bound, so deleting it deletes
    # the requirement rather than violating it, and the condition comes
    # back feasible. Take them directly.
    forced = {
        rxn_id
        for rxn_id in pool
        if constrained.reactions.get_by_id(rxn_id).lower_bound > 0.0
        or constrained.reactions.get_by_id(rxn_id).upper_bound < 0.0
    }
    pool -= forced

    previous = None
    while pool and previous != len(pool):
        previous = len(pool)
        with constrained:
            add_pfba(constrained, fraction_of_optimum=0.0)
            if not _feasible(constrained):
                break
            pool = {
                rxn_id
                for rxn_id in pool
                if abs(constrained.reactions.get_by_id(rxn_id).flux) > FLUX_EPSILON
            }

    essential = set(forced)
    for rxn_id in sorted(pool):
        with constrained:
            constrained.reactions.get_by_id(rxn_id).bounds = (0.0, 0.0)
            constrained.objective = constrained.problem.Objective(0.0)
            if not _feasible(constrained):
                essential.add(rxn_id)
    return essential


def _feasible(model: "Model") -> bool:
    """Return whether the model's constraints can be satisfied at all.

    Parameters
    ----------
    model : cobra.Model
        The model to solve, with a zero objective.

    Returns
    -------
    bool
        Whether the solver reached optimality, which with a zero objective
        means only that a feasible point exists.

    """
    model.slim_optimize()
    return model.solver.status == OPTIMAL


def _costs(
    weights: Optional[Mapping[str, float]],
    counted: Sequence[str],
    protected: Set[str],
) -> Dict[str, float]:
    """Return what each counted reaction costs to keep.

    Parameters
    ----------
    weights : dict of {str: float}, optional
        The caller's costs; anything left out costs 1.
    counted : sequence of str
        The reactions that get a binary.
    protected : set of str
        The reactions that do not, so that a weight given for one can be
        reported as having no effect rather than quietly having none.

    Returns
    -------
    dict of {str: float}
        The cost of each counted reaction.

    Raises
    ------
    SpectraError
        If any weight is negative.

    """
    costs = {rxn_id: float((weights or {}).get(rxn_id, 1.0)) for rxn_id in counted}
    negative = sorted(rxn_id for rxn_id, cost in costs.items() if cost < 0.0)
    if negative:
        raise SpectraError(
            f"weights must be non-negative, since they are what keeping a "
            f"reaction costs, but {len(negative)} are not, for example "
            f"{negative[:5]}. A negative cost would pay the solver to keep "
            f"reactions, and the smallest network is not what you would get."
        )
    ignored = sorted(set(weights or ()) - set(counted))
    if ignored:
        logger.warning(
            "%d weight(s) apply to no counted reaction and so do nothing, for "
            "example %s. A reaction in keep_reactions has no binary to put a "
            "cost on; anything else is not in the model.",
            len(ignored),
            [
                f"{rxn_id} ({'protected' if rxn_id in protected else 'absent'})"
                for rxn_id in ignored[:5]
            ],
        )
    return costs


def _check_tasks_are_possible(model: "Model", tasks: Sequence["MetabolicTask"]) -> None:
    """Raise unless the full network can already perform every task.

    Parameters
    ----------
    model : cobra.Model
        The working copy, before anything is removed.
    tasks : sequence of MetabolicTask
        The tasks required of the minimal reactome.

    Raises
    ------
    SpectraError
        If any task the full network cannot perform was required.

    Notes
    -----
    Removing reactions cannot make a task work, so a task the starting
    model fails makes the MILP infeasible -- and an infeasible MILP says
    nothing about which requirement was impossible. Checking first turns
    that into an answer. This is the mirror of
    :func:`~spectra_cobra.gapfill_for_tasks`, which skips such tasks
    because there the model is about to grow rather than shrink.

    """
    if not tasks:
        return
    from .tasks import check_tasks

    wanted = [task for task in tasks if not task.should_fail]
    impossible = [
        result.task.id for result in check_tasks(model, wanted) if not result.ok
    ]
    if impossible:
        raise SpectraError(
            f"The model cannot perform {len(impossible)} of the tasks required "
            f"of its minimal reactome, for example {impossible[:5]}. Removing "
            f"reactions cannot make a task work, so there is no solution. "
            f"Gap-fill towards them first with gapfill_for_tasks, or leave "
            f"them out."
        )


def _check_growth_is_possible(
    model: "Model",
    medium: Optional[Medium],
    biomass: Optional[str],
    required_growth: float,
    products: Mapping[str, float],
    label: str = GROWTH_CONDITION,
) -> None:
    """Raise unless the full network can meet the growth condition.

    Parameters
    ----------
    model : cobra.Model
        The working copy, before anything is removed.
    medium : dict of {str: float}, optional
        The medium.
    biomass : str, optional
        The biomass reaction.
    required_growth : float
        The rate demanded.
    products : dict of {str: float}
        The production demanded.

    Raises
    ------
    SpectraError
        If the requirements cannot all hold at once.

    Notes
    -----
    Growth and production are one flux vector, so demanding both can ask
    for more than the medium supplies even though each is reachable on
    its own. The minimal reactome is a subnetwork and can only do less,
    so the MILP comes back infeasible -- and an infeasible MILP does not
    say which requirement was the impossible one. Asking the full network
    first does.

    """
    if required_growth <= 0.0 and not products:
        return
    with _growth_condition(
        model, medium, biomass, required_growth, products
    ) as constrained:
        constrained.objective = Zero
        constrained.slim_optimize()
        feasible = constrained.solver.status == OPTIMAL
    if feasible:
        return

    wanted = []
    if required_growth > 0.0 and biomass is not None:
        wanted.append(f"{biomass} at {required_growth:.6g}")
    wanted.extend(f"{rxn_id} at {rate:.6g}" for rxn_id, rate in products.items())

    # Which of them is reachable on its own, to point at the culprit.
    singly: List[Tuple[Optional[str], float, Dict[str, float]]] = []
    if required_growth > 0.0 and biomass is not None:
        singly.append((biomass, required_growth, {}))
    singly.extend((None, 0.0, {rxn_id: rate}) for rxn_id, rate in products.items())

    alone = []
    for description, (one_biomass, one_growth, one_product) in zip(wanted, singly):
        with _growth_condition(
            model, medium, one_biomass, one_growth, one_product
        ) as constrained:
            constrained.objective = Zero
            constrained.slim_optimize()
            reachable = constrained.solver.status == OPTIMAL
        alone.append(f"{description} ({'reachable' if reachable else 'not reachable'})")

    raise SpectraError(
        f"The model cannot meet condition {label!r} even with every reaction "
        f"available, so no subnetwork of it can either. Required together in "
        f"one flux distribution: " + "; ".join(alone) + ". Growth and "
        "production share a medium, so each can be reachable on its own and "
        "the pair still not be. Lower one of them, or open the medium."
    )


def _biomass_reaction(
    model: "Model", given: Optional[str], growth_requested: bool
) -> Optional[str]:
    """Return the reaction whose flux counts as growth, if there is one.

    Parameters
    ----------
    model : cobra.Model
        The model to inspect.
    given : str, optional
        The caller's choice.
    growth_requested : bool
        Whether the caller asked for a growth rate, which makes an
        ambiguous objective an error rather than something to shrug at.

    Returns
    -------
    str or None
        The identifier, or None if there is no unambiguous one and none
        was needed.

    Raises
    ------
    SpectraError
        If the named reaction is absent, or growth was demanded of a model
        whose objective does not say which reaction it is.

    """
    if given is not None:
        if given not in model.reactions:
            raise SpectraError(f"No reaction {given!r} in the model.")
        return given
    objective = [r.id for r in model.reactions if r.objective_coefficient != 0]
    if len(objective) == 1:
        return objective[0]
    if growth_requested:
        raise SpectraError(
            f"A growth requirement was given, but the model has {len(objective)} "
            f"reactions in its objective, so which one is growth is ambiguous. "
            f"Pass biomass_reaction."
        )
    return None


def _check_medium(model: "Model", medium: Medium) -> None:
    """Raise if the medium names a reaction the model does not have.

    Parameters
    ----------
    model : cobra.Model
        The model the medium will be applied to.
    medium : dict of {str: float}
        The medium.

    Raises
    ------
    SpectraError
        If anything named is absent.

    """
    unknown = [r for r in medium if r not in model.reactions]
    if unknown:
        raise SpectraError(
            f"The medium names {len(unknown)} reaction(s) the model does not "
            f"have, for example {sorted(unknown)[:5]}."
        )


def _apply_medium(model: "Model", medium: Optional[Medium]) -> None:
    """Close every exchange, then open the ones the medium names.

    Parameters
    ----------
    model : cobra.Model
        The model to constrain, inside a ``with`` block.
    medium : dict of {str: float}, optional
        The lower bound to give each exchange to open. None leaves the
        model's own bounds alone.

    """
    if medium is None:
        return
    for reaction in model.boundary:
        reaction.lower_bound = 0.0
    for rxn_id, lower in medium.items():
        # The medium is validated against the model it was given for, but it
        # is applied again to the reduced one, which no longer has the
        # exchanges the minimisation found it did not need. An exchange that
        # is gone is closed, which is what the medium entry would have done
        # for it anyway.
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = float(lower)


def _resolve_products(
    model: "Model", products: Mapping[str, float]
) -> Dict[str, float]:
    """Turn requested products into reactions with a rate to demand.

    Parameters
    ----------
    model : cobra.Model
        The working copy. It is not modified.
    products : dict of {str: float}
        What must be made, keyed by reaction or metabolite identifier.

    Returns
    -------
    dict of {str: float}
        The flux to demand of each reaction, as a lower bound.

    Raises
    ------
    SpectraError
        If the rate is not positive, the identifier is neither a reaction
        nor a metabolite, the metabolite has no boundary reaction, it has
        more than one, or the one it has is written as a source.

    Notes
    -----
    A metabolite is resolved to the boundary reaction the model already
    has for it, and **nothing is added**. Adding a demand reaction for an
    unexchanged metabolite is the obvious convenience and it is wrong: the
    result would be the smallest subnetwork of a model the caller does not
    have, feasible only because of the drain that was added to measure it.
    Take the drain away and the same reaction set may produce nothing. So
    an unexchanged metabolite is refused, and the caller adds the reaction
    to their own model if that is really what they mean -- where it is
    visible, and counted like everything else.

    """
    resolved: Dict[str, float] = {}
    boundary = set(model.boundary)
    for key, rate in products.items():
        rate = float(rate)
        if rate <= 0.0:
            raise SpectraError(
                f"The production rate demanded of {key!r} is {rate:g}, which "
                f"requires nothing: a zero flux satisfies it. Give a positive "
                f"rate."
            )
        if key in model.reactions:
            resolved[key] = rate
            continue
        if key not in model.metabolites:
            raise SpectraError(
                f"{key!r} is neither a reaction nor a metabolite of the model, "
                f"so there is no way to require its production."
            )

        metabolite = model.metabolites.get_by_id(key)
        exchanges = sorted(rxn.id for rxn in metabolite.reactions if rxn in boundary)
        if not exchanges:
            raise SpectraError(
                f"{key!r} has no exchange reaction, so its production cannot "
                f"be required without adding one -- and adding one would "
                f"change the network being minimised. The answer would be the "
                f"smallest subnetwork of a model you do not have, feasible "
                f"only because of the drain added to measure it; without that "
                f"drain the same reactions may produce nothing at all. Name a "
                f"metabolite the model exchanges, or add the reaction to your "
                f"own model and name it."
            )
        if len(exchanges) > 1:
            raise SpectraError(
                f"{key!r} has {len(exchanges)} boundary reactions, "
                f"{exchanges}, so which one its production should be measured "
                f"through is ambiguous. Name the reaction rather than the "
                f"metabolite."
            )

        rxn_id = exchanges[0]
        coefficient = model.reactions.get_by_id(rxn_id).get_coefficient(key)
        if coefficient > 0.0:
            raise SpectraError(
                f"{rxn_id} is written as a source of {key!r} (coefficient "
                f"{coefficient:+g}), so flux through it is uptake and "
                f"requiring production would mean a negative bound rather "
                f"than a positive one. Reverse the reaction, or name a "
                f"reaction whose forward direction is the production meant."
            )
        resolved[rxn_id] = rate
        logger.info("the production of %s is required through %s", key, rxn_id)
    return resolved


def _required_growth(
    model: "Model",
    medium: Optional[Medium],
    biomass: Optional[str],
    min_growth: Optional[float],
    growth_fraction: Optional[float],
    label: str = GROWTH_CONDITION,
) -> float:
    """Return the growth rate to impose as a lower bound.

    Parameters
    ----------
    model : cobra.Model
        The working copy.
    medium : dict of {str: float}, optional
        The medium to measure the achievable rate in.
    biomass : str, optional
        The biomass reaction, if there is one.
    min_growth : float, optional
        An absolute rate.
    growth_fraction : float, optional
        A fraction of the achievable rate.
    label : str, optional
        The medium's label, for the messages.

    Returns
    -------
    float
        The rate to demand; zero means do not constrain growth.

    Raises
    ------
    SpectraError
        If a fraction was asked for of a model that cannot grow in the
        medium at all, since there is then nothing to take a fraction of.

    """
    import numpy as np

    if biomass is None:
        return 0.0
    if min_growth is not None:
        return max(float(min_growth), 0.0)

    fraction = (
        DEFAULT_GROWTH_FRACTION if growth_fraction is None else float(growth_fraction)
    )
    if fraction <= 0.0:
        return 0.0

    with model:
        _apply_medium(model, medium)
        model.objective = model.reactions.get_by_id(biomass)
        achievable = model.slim_optimize()
    if achievable is None or np.isnan(achievable) or achievable <= 0.0:
        raise SpectraError(
            f"growth_fraction={fraction:g} was asked for, but {biomass} cannot "
            f"carry flux in medium {label!r} (the model reaches "
            f"{0.0 if achievable is None or np.isnan(achievable) else achievable:g}), "
            f"so there is no rate to take a fraction of. Open the medium, or "
            f"pass min_growth instead."
        )
    logger.info(
        "%s: the model reaches %.6g; requiring %.6g (%.4g of it)",
        label,
        achievable,
        fraction * achievable,
        fraction,
    )
    return fraction * achievable


@contextmanager
def _growth_condition(
    model: "Model",
    medium: Optional[Medium],
    biomass: Optional[str],
    required_growth: float,
    products: Mapping[str, float],
) -> Iterator["Model"]:
    """Apply the medium, the growth floor and the product floors.

    Parameters
    ----------
    model : cobra.Model
        The working copy.
    medium : dict of {str: float}, optional
        The medium.
    biomass : str, optional
        The biomass reaction.
    required_growth : float
        The rate to force through it, if positive.
    products : dict of {str: float}
        The rate to force through each product reaction.

    Yields
    ------
    cobra.Model
        The same model, constrained.

    Notes
    -----
    The flux forced through `keep_reactions` is **not** applied here. It
    is the core set of the extraction, so
    :func:`~spectra_cobra.spectra_me` imposes it on the replicated model
    where the solve happens.

    """

    with model:
        _apply_medium(model, medium)
        if biomass is not None and required_growth > 0.0:
            reaction = model.reactions.get_by_id(biomass)
            reaction.lower_bound = float(required_growth)
            if reaction.upper_bound < required_growth:
                reaction.upper_bound = float(required_growth)
        for rxn_id, rate in products.items():
            reaction = model.reactions.get_by_id(rxn_id)
            reaction.lower_bound = float(rate)
            if reaction.upper_bound < rate:
                reaction.upper_bound = float(rate)
        yield model


def _can_carry(model: "Model", rxn_id: str, tol: float) -> bool:
    """Return whether a reaction can reach `tol` either way here.

    Parameters
    ----------
    model : cobra.Model
        A model already inside its condition's context.
    rxn_id : str
        The reaction to test.
    tol : float
        The flux it has to be able to reach.

    Returns
    -------
    bool
        Whether some feasible solution gives it at least `tol` in
        magnitude.

    """
    import numpy as np

    reaction = model.reactions.get_by_id(rxn_id)
    for direction in ("max", "min"):
        with model:
            model.objective = reaction
            model.objective.direction = direction
            value = model.slim_optimize()
        if value is not None and not np.isnan(value) and abs(value) >= tol:
            return True
    return False


def _forceable_keeps(
    model: "Model",
    labelled_media: Sequence[Tuple[str, Optional[Medium]]],
    biomass: Optional[str],
    required_growth: Mapping[str, float],
    products: Mapping[str, float],
    keeps: Set[str],
    tol: float,
) -> Tuple[Dict[str, Set[str]], Dict[str, Tuple[str, ...]]]:
    """Find which kept reactions can be made to carry flux where.

    Parameters
    ----------
    model : cobra.Model
        The working copy.
    labelled_media : sequence of (str, dict or None)
        The media, one growth condition each.
    biomass : str, optional
        The biomass reaction.
    required_growth : dict of {str: float}
        The rate demanded in each medium.
    products : dict of {str: float}
        What must be produced.
    keeps : set of str
        The caller's `keep_reactions`.
    tol : float
        The flux each of them has to be able to carry.

    Returns
    -------
    dict of {str: set of str}
        Per growth condition, the kept reactions that can reach `tol` in
        it and so can be made core there.
    dict of {str: tuple of str}
        Per growth condition, the ones that cannot, and so are retained
        without being forced.

    Notes
    -----
    Retaining a reaction is not the same as the reaction working. A
    reaction kept but left at zero flux is dead weight in the result, so
    `keep_reactions` becomes the **core set** of the extraction:
    :func:`~spectra_cobra.spectra_me` settles which way each one runs with
    its two LPs and then forces it to carry at least `tol`.

    Only in the **growth conditions**. A task supplies its own medium and
    closes every exchange the model has, so forcing the caller's reactions
    inside someone else's environment would be both surprising and, for
    any boundary reaction, immediately infeasible.

    This pass exists so that a kept reaction which simply cannot work in
    some medium is reported rather than turned into an infeasible core,
    which would fail the whole call and name nothing useful.

    """
    forceable: Dict[str, Set[str]] = {}
    blocked_keeps: Dict[str, Tuple[str, ...]] = {}
    if not keeps:
        return forceable, blocked_keeps

    forcing = [
        (label, one)
        for label, one in labelled_media
        if required_growth[label] > 0.0 or products
    ]
    if not forcing:
        logger.warning(
            "keep_reactions was given but there is no growth condition to "
            "force those reactions in, so they are retained without being "
            "required to carry flux. Each task defines its own medium and is "
            "left alone."
        )
        return forceable, blocked_keeps

    for label, one in forcing:
        with _growth_condition(
            model, one, biomass, required_growth[label], products
        ) as constrained:
            usable = {r for r in sorted(keeps) if _can_carry(constrained, r, tol)}
        stuck = tuple(sorted(keeps - usable))
        if usable:
            forceable[label] = usable
        if stuck:
            blocked_keeps[label] = stuck
            logger.warning(
                "%d of the %d reactions in keep_reactions cannot carry %.3g "
                "flux in condition %r, so they are retained but not forced "
                "there: %s. They will be in the result and may carry nothing.",
                len(stuck),
                len(keeps),
                tol,
                label,
                list(stuck[:10]),
            )
    return forceable, blocked_keeps


def _conditions(
    model: "Model",
    labelled_media: Sequence[Tuple[str, Optional[Medium]]],
    biomass: Optional[str],
    required_growth: Mapping[str, float],
    products: Mapping[str, float],
    tasks: Sequence["MetabolicTask"],
) -> List[Tuple[str, object]]:
    """Return the conditions to impose, each as a label and a factory.

    Parameters
    ----------
    model : cobra.Model
        The working copy, which every condition constrains in turn.
    labelled_media : sequence of (str, dict or None)
        The media, one growth condition each.
    biomass : str, optional
        The biomass reaction.
    required_growth : dict of {str: float}
        The rate to demand in each medium.
    products : dict of {str: float}
        What must be produced, in every medium.
    tasks : sequence of MetabolicTask
        The tasks, one condition each.

    Returns
    -------
    list of (str, callable)
        Each entry is a label and a no-argument callable returning a
        context manager that constrains the model for that condition.
        They are entered one at a time, so the constraints never stack.

    Raises
    ------
    SpectraError
        If two conditions would share a label.

    """
    from .tasks import task_constraints

    conditions: List[Tuple[str, object]] = []
    for label, one in labelled_media:
        # A growth condition with no floor and nothing to produce is
        # satisfied by the all-zero flux vector, so it constrains nothing
        # and only costs a copy of the network.
        if required_growth[label] > 0.0 or products:
            conditions.append(
                (
                    label,
                    lambda one=one, label=label: _growth_condition(
                        model, one, biomass, required_growth[label], products
                    ),
                )
            )
        elif one is not None:
            logger.warning(
                "medium %r requires nothing of the network -- no growth rate "
                "and no product -- so it constrains nothing and is dropped. "
                "Each task defines its own medium.",
                label,
            )

    for task in tasks:
        if task.should_fail:
            logger.warning(
                "task %r is a negative control (should_fail) and cannot be "
                "required of a minimal reactome, so it is skipped",
                task.id,
            )
            continue
        if not task.is_forced():
            logger.warning(
                "task %r forces no flux, so the all-zero solution performs it "
                "and requiring it constrains nothing",
                task.id,
            )
        conditions.append((task.id, lambda task=task: task_constraints(model, task)))

    labels = [label for label, _ in conditions]
    duplicated = sorted({label for label in labels if labels.count(label) > 1})
    if duplicated:
        raise SpectraError(
            f"Two conditions share a label, which would make them one copy of "
            f"the network rather than two: {duplicated}. Media labels and task "
            f"identifiers must all be distinct."
        )
    return conditions


def _replicate(
    base: "Model",
    conditions: Sequence[Tuple[str, object]],
    counted: Sequence[str],
) -> Tuple["Model", Dict[str, List[str]], Dict[str, Dict[str, str]]]:
    """Build one copy of the network per condition, in a single model.

    Parameters
    ----------
    base : cobra.Model
        The working copy, whose identifiers the result is reported in.
    conditions : sequence of (str, callable)
        The conditions, as returned by :func:`_conditions`.
    counted : sequence of str
        The reactions that get a binary. Anything else is free.

    Returns
    -------
    cobra.Model
        The joint model, holding every condition's copy side by side with
        no metabolite shared between them.
    dict of {str: list of str}
        For each counted reaction, the copies of it that one binary has to
        switch off together.
    dict of {str: dict of {str: str}}
        For each condition, the copy of every base reaction in it, so that
        a reaction can be named in the joint model without re-deriving the
        suffix rule anywhere else.

    Notes
    -----
    With a single condition no suffix is added and the copy is the model
    itself, so the common case pays nothing for the generality.

    """
    from cobra import Metabolite, Model, Reaction

    single = len(conditions) == 1
    joint = Model(f"{base.id}_minimal_reactome")
    joint.solver = base.problem
    joint.tolerance = base.tolerance
    groups: Dict[str, List[str]] = {rxn_id: [] for rxn_id in counted}
    copies_by_condition: Dict[str, Dict[str, str]] = {}
    # Read before any condition is entered: inside one, a task has already
    # added its own reactions to this very model, and counting those as
    # part of the network would carry the scaffolding into the answer.
    base_ids = {rxn.id for rxn in base.reactions}

    for label, factory in conditions:
        suffix = "" if single else f"{CONDITION_SEPARATOR}{label}"
        with factory() as constrained:
            metabolites = {
                met.id: Metabolite(
                    f"{met.id}{suffix}",
                    formula=met.formula,
                    name=met.name,
                    compartment=met.compartment,
                )
                for met in constrained.metabolites
            }
            joint.add_metabolites(list(metabolites.values()))
            copies = [
                Reaction(
                    f"{rxn.id}{suffix}",
                    name=rxn.name,
                    lower_bound=rxn.lower_bound,
                    upper_bound=rxn.upper_bound,
                )
                for rxn in constrained.reactions
            ]
            joint.add_reactions(copies)
            copies_by_condition[label] = {
                rxn.id: copy.id
                for rxn, copy in zip(constrained.reactions, copies)
                if rxn.id in base_ids
            }
            for rxn, copy in zip(constrained.reactions, copies):
                copy.add_metabolites(
                    {
                        metabolites[met.id]: coeff
                        for met, coeff in rxn.metabolites.items()
                    }
                )
                if rxn.id in groups:
                    groups[rxn.id].append(copy.id)

    # Every condition is built from the same model and none of them removes
    # a reaction, so today no reaction ends up without a copy. If one ever
    # did, a binary governing nothing would be free to sit at one and inflate
    # the count, so it is dropped instead -- a reaction no condition contains
    # cannot be needed by any of them.
    empty = sorted(rxn_id for rxn_id, copies in groups.items() if not copies)
    for rxn_id in empty:
        del groups[rxn_id]
    if empty:
        logger.warning(
            "%d reaction(s) appear in no condition, so nothing requires them "
            "and they are dropped, for example %s",
            len(empty),
            empty[:5],
        )
    return joint, groups, copies_by_condition


def _report(
    base: "Model",
    original_ids: Set[str],
    kept: Set[str],
    labelled_media: Sequence[Tuple[str, Optional[Medium]]],
    biomass: Optional[str],
    required_growth: Mapping[str, float],
    required_products: Mapping[str, float],
    tasks: Sequence["MetabolicTask"],
    conditions: Tuple[str, ...],
    essential: Tuple[str, ...],
    forceable: Mapping[str, Set[str]],
    kept_but_blocked: Mapping[str, Tuple[str, ...]],
    tol: float,
    check: bool,
) -> MinimalReactome:
    """Cut the model down to the kept reactions and measure what it does.

    Parameters
    ----------
    base : cobra.Model
        The working copy to cut down, already without any blocked
        reactions preprocessing deleted.
    original_ids : set of str
        Every reaction the model started with, so that what preprocessing
        deleted is still reported as removed.
    kept : set of str
        The reactions to retain.
    labelled_media : sequence of (str, dict or None)
        The media, for the verification.
    biomass : str, optional
        The biomass reaction.
    required_growth : dict of {str: float}
        The rate demanded in each medium.
    required_products : dict of {str: float}
        The production demanded.
    tasks : sequence of MetabolicTask
        The tasks that were demanded.
    conditions : tuple of str
        The condition labels, for the record.
    essential : tuple of str
        What preprocessing proved is in every answer.
    forceable : dict of {str: set of str}
        Per growth condition, the kept reactions the solve was told to
        force, so that the verification can check they really work.
    kept_but_blocked : dict of {str: tuple of str}
        Per growth condition, the kept reactions nothing could force.
    tol : float
        The flux a kept reaction has to carry.
    check : bool
        Whether to verify at all.

    Returns
    -------
    MinimalReactome
        The result.

    """
    # Metabolites left without a reaction are kept. They cost nothing, and
    # dropping them would make a task that names one impossible to even set
    # up against the result -- reported as a failure of a network that in
    # fact performs the task, since the task brings its own exchanges.
    model = subset_model(base, set(kept), keep_orphan_metabolites=True)
    removed = tuple(sorted(original_ids - kept))

    growth: Dict[str, float] = {}
    production: Dict[str, Dict[str, float]] = {}
    unsatisfied: List[str] = []
    if check:
        for label, one in labelled_media:
            if label not in conditions:
                continue
            reached, made, ok = _measure(
                model, one, biomass, required_growth[label], required_products
            )
            if reached is not None:
                growth[label] = reached
            production[label] = made
            if not ok:
                unsatisfied.append(label)
                continue
            # A forced reaction that cannot carry flux in the result would
            # mean the extraction did not keep the promise it was given.
            dead = _keeps_without_flux(
                model,
                one,
                biomass,
                required_growth[label],
                required_products,
                forceable.get(label, set()),
                tol,
            )
            if dead:
                logger.warning(
                    "%d reaction(s) forced to carry flux in condition %r "
                    "cannot do so in the result: %s",
                    len(dead),
                    label,
                    sorted(dead)[:5],
                )
                unsatisfied.append(label)
        if tasks:
            from .tasks import check_tasks

            wanted = [task for task in tasks if not task.should_fail]
            unsatisfied.extend(
                result.task.id for result in check_tasks(model, wanted) if not result.ok
            )
        if unsatisfied:
            logger.warning(
                "the minimal reactome does not meet %d of its requirements: %s. "
                "This should not happen; it usually means the solver stopped at "
                "its time limit, or that model.tolerance is loose enough for a "
                "requirement to be met only to within it.",
                len(unsatisfied),
                unsatisfied,
            )

    return MinimalReactome(
        model=model,
        kept=tuple(sorted(kept)),
        removed=removed,
        growth=growth,
        required_growth={
            label: rate
            for label, rate in required_growth.items()
            if label in conditions
        },
        production=production,
        required_production=dict(required_products),
        conditions=conditions,
        essential=essential,
        kept_but_blocked=dict(kept_but_blocked),
        unsatisfied=tuple(unsatisfied),
    )


def _measure(
    model: "Model",
    medium: Optional[Medium],
    biomass: Optional[str],
    required_growth: float,
    required_products: Mapping[str, float],
) -> Tuple[Optional[float], Dict[str, float], bool]:
    """Solve one growth condition and read off what it achieves.

    Parameters
    ----------
    model : cobra.Model
        The minimal reactome.
    medium : dict of {str: float}, optional
        The medium.
    biomass : str, optional
        The biomass reaction.
    required_growth : float
        The rate demanded.
    required_products : dict of {str: float}
        The production demanded.

    Returns
    -------
    float or None
        The growth rate reached, or None if there is no biomass reaction
        or the condition came back infeasible.
    dict of {str: float}
        The flux each required product carries in that same solution.
    bool
        Whether every requirement was met.

    Notes
    -----
    One solve, not one per requirement, because the requirements were
    imposed together and have to be met together: a network that can make
    the product *or* grow, but not both at once, has not met them.

    """
    import numpy as np

    present = {rxn.id for rxn in model.reactions}
    missing = sorted(set(required_products) - present)
    if missing or (biomass is not None and biomass not in present):
        return None, {}, False

    with _growth_condition(
        model, medium, biomass, required_growth, required_products
    ) as constrained:
        if biomass is not None:
            constrained.objective = constrained.reactions.get_by_id(biomass)
        value = constrained.slim_optimize()
        if value is None or np.isnan(value):
            return None, {}, False
        growth = (
            float(constrained.reactions.get_by_id(biomass).flux)
            if biomass is not None
            else None
        )
        production = {
            rxn_id: float(constrained.reactions.get_by_id(rxn_id).flux)
            for rxn_id in required_products
        }

    floor = 1.0 - VERIFICATION_SLACK
    ok = growth is None or growth >= required_growth * floor
    for rxn_id, rate in required_products.items():
        ok = ok and production[rxn_id] >= rate * floor
    return growth, production, ok


def _keeps_without_flux(
    model: "Model",
    medium: Optional[Medium],
    biomass: Optional[str],
    required_growth: float,
    required_products: Mapping[str, float],
    forced: Set[str],
    tol: float,
) -> Set[str]:
    """Return the forced reactions the result cannot make work.

    Parameters
    ----------
    model : cobra.Model
        The minimal reactome.
    medium : dict of {str: float}, optional
        The medium of this growth condition.
    biomass : str, optional
        The biomass reaction.
    required_growth : float
        The rate demanded.
    required_products : dict of {str: float}
        The production demanded.
    forced : set of str
        The kept reactions the solve was told to force here.
    tol : float
        The flux each of them has to be able to carry.

    Returns
    -------
    set of str
        The ones that cannot. Expected to be empty: the extraction forced
        them, so the result must be able to.

    Notes
    -----
    Checked by outcome rather than by re-imposing the bound. Re-imposing
    it would only confirm that the same constraint is still satisfiable,
    which it is by construction; asking whether the reaction can carry
    flux in the network that came back is the question worth asking.

    """
    if not forced:
        return set()
    with _growth_condition(
        model, medium, biomass, required_growth, required_products
    ) as constrained:
        return {
            rxn_id
            for rxn_id in sorted(forced)
            if rxn_id in constrained.reactions
            and not _can_carry(constrained, rxn_id, tol)
        }
