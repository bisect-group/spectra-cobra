"""Metabolic tasks, and the reactions that are essential to them.

A *task* is a statement about what a network must be able to do: take up
these metabolites, produce those, and carry flux through this reaction. It is
checked by closing the model's boundary, opening exchanges for the task's
inputs and outputs, adding any reactions the task defines, and asking whether
the result is feasible.

Tasks matter for extraction because of one observation. Requiring a model to
perform a task is awkward to express as a constraint, but requiring it to
*contain a reaction* is exactly what :func:`spectra_cobra.spectra_me` already
does. So compute, for each task, the reactions without which the task becomes
infeasible, and hand their union to the extraction as core reactions.

That is a necessary condition rather than a sufficient one. A task with two
alternative routes has no essential reactions at all, since either route can
be removed on its own, so forcing the essential set does not guarantee the
task survives. It carries most of the weight cheaply; :func:`check_tasks` on
the extracted model is what tells you whether anything is left to repair.

The approach is taken from ftINIT, which computes the reactions essential to
each task and forces them to carry flux in the extraction that follows, then
gap-fills task by task for whatever is left over:

    Gustafsson, J., Anton, M., Roshanzamir, F., Jörnsten, R., Kerkhoven,
    E. J., Robinson, J. L., and Nielsen, J. (2023). Generation and analysis
    of context-specific genome-scale metabolic models derived from
    single-cell RNA-Seq data. *Proceedings of the National Academy of
    Sciences*, 120(6), e2217868120. https://doi.org/10.1073/pnas.2217868120

The task file format :func:`parse_task_list` reads, and the published task
lists written in it, come from the same lineage by way of the RAVEN Toolbox.
"""

import re
from contextlib import contextmanager
from dataclasses import dataclass, field
from logging import getLogger
from typing import TYPE_CHECKING, Dict, Iterable, Iterator, List, Optional, Set, Tuple

from optlang.interface import OPTIMAL

from .exceptions import SpectraError

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)

#: The bounds a task gives an input or output exchange when it does not say,
#: following the convention of the published task lists: uptake and production
#: are allowed up to 1000 but never required.
DEFAULT_BOUNDS = (0.0, 1000.0)

#: The bounds a task equation gets when it does not say.
DEFAULT_EQUATION_BOUNDS = (-1000.0, 1000.0)

#: Below this absolute flux a reaction is treated as unused when picking the
#: candidates for the essentiality search.
FLUX_EPSILON = 1e-9


@dataclass(frozen=True)
class TaskEquation:
    """A reaction a task adds to the model while it is being checked.

    Parameters
    ----------
    reaction : dict of {str: float}
        The stoichiometry, keyed by metabolite identifier. Negative
        coefficients consume, positive produce.
    bounds : tuple of (float, float), optional
        The flux bounds (default ``(-1000, 1000)``). A positive lower bound is
        what forces the task to actually happen.

    """

    reaction: Dict[str, float]
    bounds: Tuple[float, float] = DEFAULT_EQUATION_BOUNDS


@dataclass(frozen=True)
class MetabolicTask:
    """Something a metabolic network must be able to do.

    Parameters
    ----------
    id : str
        A short identifier.
    description : str, optional
        What the task is for, in words.
    inputs : dict of {str: tuple of (float, float)}, optional
        The metabolites that may be taken up, keyed by identifier, each with
        the bounds on its uptake flux (default ``(0, 1000)`` each, i.e.
        allowed but not required).
    outputs : dict of {str: tuple of (float, float)}, optional
        The metabolites that may be produced, with the bounds on their
        production flux. A positive lower bound *requires* the production.
    equations : iterable of TaskEquation, optional
        Reactions added for the duration of the check.
    should_fail : bool, optional
        Whether the task is expected to be infeasible (default False). These
        are the negative controls of a task list: a model that can do them is
        leaking.

    Notes
    -----
    A task with no positive lower bound anywhere is satisfied by the
    all-zero flux vector, so it passes trivially. Force the task through an
    output's lower bound or a :class:`TaskEquation`'s.

    """

    id: str
    description: str = ""
    inputs: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    outputs: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    equations: Tuple[TaskEquation, ...] = ()
    should_fail: bool = False

    def is_forced(self) -> bool:
        """Return whether anything in the task requires nonzero flux.

        Returns
        -------
        bool
            Whether some output or equation has a positive lower bound.

        """
        forced_out = any(lower > 0.0 for lower, _ in self.outputs.values())
        forced_equ = any(equation.bounds[0] > 0.0 for equation in self.equations)
        return forced_out or forced_equ


@dataclass(frozen=True)
class TaskResult:
    """The outcome of checking one task.

    Parameters
    ----------
    task : MetabolicTask
        The task that was checked.
    feasible : bool
        Whether the constrained model had a solution.
    ok : bool
        Whether that is the expected outcome, i.e. `feasible` unless the task
        was marked ``should_fail``.
    missing : tuple of str
        Identifiers the task named that the model does not have. A task
        naming a missing metabolite is reported rather than checked.

    """

    task: MetabolicTask
    feasible: bool
    ok: bool
    missing: Tuple[str, ...] = ()


def _resolve(model: "Model", identifiers: Iterable[str]) -> List[str]:
    """Return the identifiers the model does not have.

    Parameters
    ----------
    model : cobra.Model
        The model to look in.
    identifiers : iterable of str
        Metabolite identifiers.

    Returns
    -------
    list of str
        Those that are absent.

    """
    present = {met.id for met in model.metabolites}
    return [met_id for met_id in identifiers if met_id not in present]


def missing_metabolites(model: "Model", task: "MetabolicTask") -> List[str]:
    """Return the metabolites a task names that a model does not have.

    Parameters
    ----------
    model : cobra.Model
        The model to check against.
    task : MetabolicTask
        The task to inspect.

    Returns
    -------
    list of str
        The absent metabolite identifiers, in the order inputs, outputs,
        equations.

    """
    named = list(task.inputs) + list(task.outputs)
    for equation in task.equations:
        named.extend(equation.reaction)
    seen, ordered = set(), []
    for met_id in named:
        if met_id not in seen:
            seen.add(met_id)
            ordered.append(met_id)
    return _resolve(model, ordered)


@contextmanager
def task_constraints(model: "Model", task: "MetabolicTask") -> Iterator["Model"]:
    """Apply a task to a model for the duration of the block.

    Every existing boundary reaction is closed, exchanges are opened for the
    task's inputs and outputs, and the task's equations are added. Everything
    is undone on the way out.

    Parameters
    ----------
    model : cobra.Model
        The model to constrain.
    task : MetabolicTask
        The task to apply.

    Yields
    ------
    cobra.Model
        The same model, constrained.

    Raises
    ------
    SpectraError
        If the task names a metabolite the model does not have.

    """
    from cobra import Reaction

    absent = missing_metabolites(model, task)
    if absent:
        raise SpectraError(
            f"Task {task.id!r} names {len(absent)} metabolite(s) the model does "
            f"not have, for example {absent[:5]}."
        )

    with model:
        # The task defines its own medium, so the model's own exchanges must
        # not be able to supply anything alongside it.
        for reaction in model.boundary:
            reaction.bounds = (0.0, 0.0)

        added = []
        for index, (met_id, (lower, upper)) in enumerate(task.inputs.items()):
            # A supply reaction, written as "--> met", so its flux is the
            # uptake rate directly rather than the negated rate an exchange
            # written "met -->" would carry.
            reaction = Reaction(f"spectra_task_in_{index}_{met_id}")
            reaction.bounds = (float(lower), float(upper))
            added.append((reaction, {met_id: 1.0}))
        for index, (met_id, (lower, upper)) in enumerate(task.outputs.items()):
            reaction = Reaction(f"spectra_task_out_{index}_{met_id}")
            reaction.bounds = (float(lower), float(upper))
            added.append((reaction, {met_id: -1.0}))
        for index, equation in enumerate(task.equations):
            reaction = Reaction(f"spectra_task_equ_{index}")
            reaction.bounds = (float(equation.bounds[0]), float(equation.bounds[1]))
            added.append((reaction, dict(equation.reaction)))

        model.add_reactions([reaction for reaction, _ in added])
        for reaction, stoichiometry in added:
            reaction.add_metabolites(
                {model.metabolites.get_by_id(m): c for m, c in stoichiometry.items()}
            )
        yield model


def _is_feasible(model: "Model") -> bool:
    """Solve for feasibility alone and report whether a solution exists.

    Parameters
    ----------
    model : cobra.Model
        The model to solve, with its objective already zeroed.

    Returns
    -------
    bool
        Whether the solver reached optimality, which with a zero objective
        means nothing more than that the constraints can be satisfied.

    Notes
    -----
    ``slim_optimize(error_value=None)`` raises on a non-optimal status rather
    than returning None, which is the wrong shape when infeasibility is the
    answer being measured, so the status is read directly instead.

    """
    model.slim_optimize()
    return model.solver.status == OPTIMAL


def check_task(model: "Model", task: "MetabolicTask") -> TaskResult:
    """Check whether a model can perform a task.

    Parameters
    ----------
    model : cobra.Model
        The model to check.
    task : MetabolicTask
        The task to check.

    Returns
    -------
    TaskResult
        The outcome. A task naming metabolites the model lacks is reported as
        infeasible with those identifiers in ``missing`` rather than raising.

    """
    absent = missing_metabolites(model, task)
    if absent:
        return TaskResult(task, False, task.should_fail, tuple(absent))

    with task_constraints(model, task):
        model.objective = model.problem.Objective(0.0)
        feasible = _is_feasible(model)
    return TaskResult(task, feasible, feasible != task.should_fail)


def check_tasks(model: "Model", tasks: Iterable["MetabolicTask"]) -> List[TaskResult]:
    """Check a model against a list of tasks.

    Parameters
    ----------
    model : cobra.Model
        The model to check.
    tasks : iterable of MetabolicTask
        The tasks to check.

    Returns
    -------
    list of TaskResult
        One result per task, in the order given.

    """
    results = [check_task(model, task) for task in tasks]
    failed = [r.task.id for r in results if not r.ok]
    logger.info(
        "%d of %d tasks behaved as expected%s",
        len(results) - len(failed),
        len(results),
        "" if not failed else f"; unexpected: {failed[:5]}",
    )
    return results


def essential_reactions_for_task(
    model: "Model",
    task: "MetabolicTask",
    candidates: Optional[Iterable[str]] = None,
) -> Set[str]:
    """Find the reactions without which a task becomes infeasible.

    Parameters
    ----------
    model : cobra.Model
        The model to search in.
    task : MetabolicTask
        The task in question. A task the model cannot perform at all has no
        essential reactions.
    candidates : iterable of str, optional
        Restrict the search to these reactions (default: all of them).

    Returns
    -------
    set of str
        The identifiers of the essential reactions.

    Notes
    -----
    Only a reaction carrying flux in a feasible solution can be essential: if
    some solution leaves it at zero, fixing it to zero keeps that solution
    available. The search therefore starts from the reactions used by one
    solution and shrinks that set by repeatedly minimising the flux through
    whatever is left, which drops in one LP every reaction that some other
    route can avoid. Only the survivors are then tested one at a time.

    """
    from cobra.flux_analysis.parsimonious import add_pfba

    absent = missing_metabolites(model, task)
    if absent:
        return set()

    with task_constraints(model, task) as constrained:
        pool = (
            {rxn.id for rxn in constrained.reactions}
            if candidates is None
            else set(candidates)
        )
        # Never report the scaffolding this function added as essential.
        pool = {r for r in pool if not r.startswith("spectra_task_")}

        constrained.objective = constrained.problem.Objective(0.0)
        if not _is_feasible(constrained):
            return set()

        # Shrink the candidates: anything a minimal-flux solution can leave
        # at zero is not essential, and one LP rules out many at once.
        previous = None
        while previous != len(pool):
            previous = len(pool)
            with constrained:
                add_pfba(constrained, fraction_of_optimum=0.0)
                if not _is_feasible(constrained):
                    break
                pool = {
                    rxn_id
                    for rxn_id in pool
                    if abs(constrained.reactions.get_by_id(rxn_id).flux) > FLUX_EPSILON
                }

        essential = set()
        for rxn_id in sorted(pool):
            reaction = constrained.reactions.get_by_id(rxn_id)
            with constrained:
                reaction.bounds = (0.0, 0.0)
                constrained.objective = constrained.problem.Objective(0.0)
                if not _is_feasible(constrained):
                    essential.add(rxn_id)
    return essential


def essential_reactions_for_tasks(
    model: "Model",
    tasks: Iterable["MetabolicTask"],
    candidates: Optional[Iterable[str]] = None,
) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """Find the reactions essential to any of a list of tasks.

    This is the step that turns "the model must perform these tasks" into a
    core reaction set that :func:`spectra_cobra.spectra_me` can enforce.

    Parameters
    ----------
    model : cobra.Model
        The model to search in, normally the universal model being extracted
        from.
    tasks : iterable of MetabolicTask
        The tasks. Those marked ``should_fail`` are skipped, since nothing
        should be forced in to support them.
    candidates : iterable of str, optional
        Restrict the search to these reactions (default: all of them).

    Returns
    -------
    tuple of (set of str, dict of {str: set of str})
        The union of the essential reactions, and the essential set of each
        task keyed by task identifier.

    """
    per_task: Dict[str, Set[str]] = {}
    for task in tasks:
        if task.should_fail:
            continue
        per_task[task.id] = essential_reactions_for_task(model, task, candidates)
        logger.debug("task %s: %d essential reactions", task.id, len(per_task[task.id]))
    union: Set[str] = set().union(*per_task.values()) if per_task else set()
    logger.info(
        "%d reactions are essential to at least one of %d tasks",
        len(union),
        len(per_task),
    )
    return union, per_task


#: How a published task file refers to a metabolite: its *name*, then its
#: compartment in brackets, as in ``glucose[e]``.
_REFERENCE = re.compile(r"^\s*(?P<name>.+?)\s*\[\s*(?P<compartment>[^\]]+)\s*\]\s*$")

#: The arrows a task equation may be written with.
_ARROWS = ("<=>", "=>", "-->", "<==>", "->")


def _metabolite_index(model: "Model") -> Dict[Tuple[str, str], str]:
    """Index a model's metabolites by (name, compartment).

    Parameters
    ----------
    model : cobra.Model
        The model to index.

    Returns
    -------
    dict of {(str, str): str}
        The identifier of each metabolite, keyed by its name and
        compartment. Where several share both, the first is kept.

    """
    index: Dict[Tuple[str, str], str] = {}
    for metabolite in model.metabolites:
        key = (metabolite.name, metabolite.compartment)
        index.setdefault(key, metabolite.id)
        # Published lists are not always careful about case.
        index.setdefault(
            (metabolite.name.lower(), metabolite.compartment), metabolite.id
        )
    return index


def _resolve_reference(
    reference: str, index: Dict[Tuple[str, str], str]
) -> Optional[str]:
    """Return the identifier a ``name[compartment]`` reference points to.

    Parameters
    ----------
    reference : str
        The reference, as written in the task file.
    index : dict of {(str, str): str}
        The index from :func:`_metabolite_index`.

    Returns
    -------
    str or None
        The metabolite identifier, or None if it is not in the model.

    """
    match = _REFERENCE.match(reference)
    if match is None:
        return None
    name = match.group("name")
    compartment = match.group("compartment")
    return index.get((name, compartment)) or index.get((name.lower(), compartment))


def _parse_equation(
    equation: str, index: Dict[Tuple[str, str], str]
) -> Tuple[Dict[str, float], List[str]]:
    """Turn a task equation into stoichiometry.

    Parameters
    ----------
    equation : str
        The equation, such as ``ATP[c] + H2O[c] => ADP[c] + Pi[c]``.
    index : dict of {(str, str): str}
        The metabolite index.

    Returns
    -------
    tuple of (dict of {str: float}, list of str)
        The stoichiometry, and the references that could not be resolved.

    """
    arrow = next((a for a in _ARROWS if a in equation), None)
    if arrow is None:
        return {}, [equation]
    left, right = equation.split(arrow, 1)
    stoichiometry: Dict[str, float] = {}
    unresolved: List[str] = []
    for side, sign in ((left, -1.0), (right, 1.0)):
        # Split on " + ", not on "+": the published lists write the proton
        # as "H+[c]", and a bare split shears it into "H" and "[c]".
        for term in re.split(r"\s\+\s", side):
            term = term.strip()
            if not term:
                continue
            coefficient, _, remainder = term.partition(" ")
            try:
                factor = float(coefficient)
                reference = remainder.strip()
            except ValueError:
                factor = 1.0
                reference = term
            met_id = _resolve_reference(reference, index)
            if met_id is None:
                unresolved.append(reference)
                continue
            stoichiometry[met_id] = stoichiometry.get(met_id, 0.0) + sign * factor
    return stoichiometry, unresolved


def _bounds(row: Dict[str, str], prefix: str, default: Tuple[float, float]):
    """Read a row's lower and upper bound for a column group.

    Parameters
    ----------
    row : dict of {str: str}
        The task row.
    prefix : str
        The column prefix, ``"IN"``, ``"OUT"`` or ``"EQU"``.
    default : tuple of (float, float)
        What to use where the row is blank.

    Returns
    -------
    tuple of (float, float)
        The bounds.

    """
    lower, upper = default
    raw_lower = (row.get(f"{prefix} LB") or "").strip()
    raw_upper = (row.get(f"{prefix} UB") or "").strip()
    if raw_lower:
        lower = float(raw_lower)
    if raw_upper:
        upper = float(raw_upper)
    return lower, upper


def parse_task_list(
    path: str,
    model: "Model",
    skip_unresolved: bool = True,
    free_outputs: Iterable[str] = (),
) -> Tuple[List[MetabolicTask], Dict[str, List[str]]]:
    """Read a published task list and bind it to a model.

    Parameters
    ----------
    path : str
        The tab-separated task file. The columns follow the convention of
        the published lists: ``ID``, ``DESCRIPTION``, ``SHOULD FAIL``,
        ``IN``, ``IN LB``, ``IN UB``, ``OUT``, ``OUT LB``, ``OUT UB``,
        ``EQU``, ``EQU LB``, ``EQU UB``, with metabolites written as
        ``name[compartment]`` and several separated by semicolons.
    model : cobra.Model
        The model the names are resolved against, by metabolite name and
        compartment rather than by identifier.
    skip_unresolved : bool, optional
        Whether to drop tasks naming metabolites the model does not have
        (default True). They can never be performed, so keeping them only
        produces failures that no gap-fill could repair.
    free_outputs : iterable of str, optional
        References, written as ``name[compartment]``, that every task may
        produce whether or not it lists them (default none). Published
        lists are often written against a convention where common waste can
        always leave, and omit it; see the notes.

    Returns
    -------
    tuple of (list of MetabolicTask, dict of {str: list of str})
        The tasks, and the references that could not be resolved, keyed by
        task identifier.

    Notes
    -----
    Published lists refer to metabolites by name, not by identifier, which
    is what makes one list usable across models that number their
    metabolites differently. It also makes the binding approximate: a name
    that does not appear verbatim in the model will not resolve, and the
    second return value is there so that is visible rather than silent.

    `free_outputs` exists because the lists are not self-contained. A task
    reading "produce 3-phospho-D-glycerate from glucose, oxygen and
    phosphate" lists the product and nothing else, but making it also makes
    water and protons, and with nowhere for those to go the task is
    infeasible as written. On Human-GEM, 30 of 56 published essential tasks
    fail for exactly that reason and pass once ``H2O``, ``CO2`` and ``H+``
    are allowed out. Which metabolites a list assumes is a property of the
    list, so it is asked for rather than guessed.

    """
    import csv

    index = _metabolite_index(model)
    tasks: List[MetabolicTask] = []
    unresolved: Dict[str, List[str]] = {}
    seen: Dict[str, int] = {}

    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            identifier = (row.get("ID") or "").strip()
            description = (row.get("DESCRIPTION") or "").strip()
            if not identifier or identifier.startswith("#") or not description:
                continue
            # The published lists reuse a category as the ID, so make it unique.
            seen[identifier] = seen.get(identifier, 0) + 1
            task_id = f"{identifier}_{seen[identifier]:03d}"

            missing: List[str] = []
            inputs, outputs = {}, {}
            for column, target, default in (
                ("IN", inputs, DEFAULT_BOUNDS),
                ("OUT", outputs, DEFAULT_BOUNDS),
            ):
                bounds = _bounds(row, column, default)
                for reference in (row.get(column) or "").split(";"):
                    reference = reference.strip()
                    if not reference:
                        continue
                    met_id = _resolve_reference(reference, index)
                    if met_id is None:
                        missing.append(reference)
                        continue
                    target[met_id] = bounds

            for reference in free_outputs:
                met_id = _resolve_reference(reference, index)
                if met_id is None:
                    missing.append(reference)
                elif met_id not in outputs:
                    outputs[met_id] = DEFAULT_BOUNDS

            equations = []
            equ_bounds = _bounds(row, "EQU", DEFAULT_EQUATION_BOUNDS)
            for raw in (row.get("EQU") or "").split(";"):
                raw = raw.strip()
                if not raw:
                    continue
                stoichiometry, bad = _parse_equation(raw, index)
                missing.extend(bad)
                if stoichiometry:
                    equations.append(TaskEquation(stoichiometry, equ_bounds))

            if missing:
                unresolved[task_id] = missing
                if skip_unresolved:
                    continue

            tasks.append(
                MetabolicTask(
                    id=task_id,
                    description=description,
                    inputs=inputs,
                    outputs=outputs,
                    equations=tuple(equations),
                    should_fail=bool((row.get("SHOULD FAIL") or "").strip()),
                )
            )

    logger.info(
        "parsed %d tasks; %d named metabolites the model does not have",
        len(tasks),
        len(unresolved),
    )
    return tasks, unresolved
