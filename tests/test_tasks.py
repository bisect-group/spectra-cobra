"""Test metabolic tasks and the essential reactions behind them."""

import pytest
from cobra import Model

from spectra_cobra import (
    MetabolicTask,
    SpectraError,
    TaskEquation,
    check_task,
    check_tasks,
    essential_reactions_for_task,
    essential_reactions_for_tasks,
)

#: Hydrolysing ATP, forced, which is how the published task lists make a
#: model demonstrate that it can actually regenerate ATP.
ATP_HYDROLYSIS = TaskEquation(
    {"atp_c": -1.0, "h2o_c": -1.0, "adp_c": 1.0, "pi_c": 1.0, "h_c": 1.0},
    (1.0, 1000.0),
)


@pytest.fixture(scope="function")
def task_model(solver: str) -> Model:
    """Return a toy model that can make ``atp_c`` from ``s_e``, one way only.

    ``T1`` imports the substrate, ``R1`` phosphorylates ADP at the cost of
    consuming it, and ``T2`` exports the waste. Nothing is redundant, so
    every reaction is essential to any task that needs ATP.

    Returns
    -------
    cobra.Model
        The toy model.

    """
    from cobra import Reaction

    model = Model("task toy")
    specs = [
        ("T1", "--> s_e", 0.0, 1000.0),
        ("T2", "s_c --> w_c", 0.0, 1000.0),
        ("T3", "w_c -->", 0.0, 1000.0),
        ("UP", "s_e --> s_c", 0.0, 1000.0),
        ("R1", "s_c + adp_c + pi_c + h_c --> atp_c + h2o_c + w_c", 0.0, 1000.0),
    ]
    model.add_reactions([Reaction(i) for i, _, _, _ in specs])
    for rxn_id, equation, lower, upper in specs:
        reaction = model.reactions.get_by_id(rxn_id)
        reaction.reaction = equation
        reaction.bounds = (lower, upper)
    model.solver = solver
    return model


def _atp_task(**kwargs) -> MetabolicTask:
    """Return the "make ATP from the substrate" task.

    Returns
    -------
    MetabolicTask
        The task.

    """
    return MetabolicTask(
        id="atp",
        description="Regenerate ATP from the substrate",
        inputs={
            "s_e": (0.0, 1000.0),
            "adp_c": (0.0, 1000.0),
            "pi_c": (0.0, 1000.0),
            "h_c": (0.0, 1000.0),
        },
        outputs={
            "w_c": (0.0, 1000.0),
            "h2o_c": (0.0, 1000.0),
            # The hydrolysis equation releases a proton, which has to go
            # somewhere or the task is infeasible at steady state.
            "h_c": (0.0, 1000.0),
        },
        equations=(ATP_HYDROLYSIS,),
        **kwargs,
    )


def test_a_task_the_model_can_do_passes(task_model: Model) -> None:
    """The toy model can regenerate ATP from its substrate."""
    result = check_task(task_model, _atp_task())

    assert result.feasible
    assert result.ok


def test_a_task_needs_forcing_to_mean_anything(task_model: Model) -> None:
    """Without a positive lower bound a task is met by doing nothing."""
    idle = MetabolicTask(id="idle", inputs={"s_e": (0.0, 1000.0)})

    assert not idle.is_forced()
    assert _atp_task().is_forced()
    # It still "passes", which is exactly why is_forced is worth checking.
    assert check_task(task_model, idle).feasible


def test_a_task_without_its_input_fails(task_model: Model) -> None:
    """Cutting off the substrate makes the ATP task infeasible."""
    starved = MetabolicTask(
        id="atp_no_substrate",
        inputs={
            "adp_c": (0.0, 1000.0),
            "pi_c": (0.0, 1000.0),
            "h_c": (0.0, 1000.0),
        },
        outputs={
            "w_c": (0.0, 1000.0),
            "h2o_c": (0.0, 1000.0),
            "h_c": (0.0, 1000.0),
        },
        equations=(ATP_HYDROLYSIS,),
    )

    assert not check_task(task_model, starved).feasible


def test_should_fail_inverts_the_expectation(task_model: Model) -> None:
    """A negative control is ok precisely when it is infeasible."""
    control = MetabolicTask(
        id="free_atp",
        equations=(ATP_HYDROLYSIS,),
        should_fail=True,
    )

    result = check_task(task_model, control)

    assert not result.feasible
    assert result.ok


def test_the_model_is_left_alone(task_model: Model) -> None:
    """Checking a task must not leave its scaffolding behind."""
    before = {rxn.id for rxn in task_model.reactions}
    bounds = {rxn.id: rxn.bounds for rxn in task_model.reactions}

    check_task(task_model, _atp_task())

    assert {rxn.id for rxn in task_model.reactions} == before
    assert {rxn.id: rxn.bounds for rxn in task_model.reactions} == bounds


def test_a_task_naming_an_unknown_metabolite_is_reported(task_model: Model) -> None:
    """Naming a metabolite the model lacks is a result, not a crash."""
    result = check_task(
        task_model, MetabolicTask(id="bogus", inputs={"nope_c": (0.0, 1.0)})
    )

    assert not result.feasible
    assert result.missing == ("nope_c",)


def test_task_constraints_raise_on_an_unknown_metabolite(task_model: Model) -> None:
    """Applying such a task directly is an error, since it cannot be built."""
    from spectra_cobra.tasks import task_constraints

    with pytest.raises(SpectraError, match="does not have"):
        with task_constraints(
            task_model, MetabolicTask("x", inputs={"nope_c": (0, 1)})
        ):
            pass  # pragma: no cover


def test_essential_reactions_are_the_ones_the_task_cannot_lose(
    task_model: Model,
) -> None:
    """Every reaction on the only route to ATP is essential to the task."""
    essential = essential_reactions_for_task(task_model, _atp_task())

    assert essential == {"UP", "R1"}


def test_a_redundant_route_leaves_nothing_essential(task_model: Model) -> None:
    """Two ways to do the same thing means neither is essential on its own.

    This is the limit of the approach: forcing the essential reactions in
    cannot guarantee a task survives, because a task with alternatives has
    no essential reactions to force.
    """
    from cobra import Reaction

    bypass = Reaction("UP2")
    task_model.add_reactions([bypass])
    bypass.reaction = "s_e --> s_c"
    bypass.bounds = (0.0, 1000.0)

    essential = essential_reactions_for_task(task_model, _atp_task())

    assert "UP" not in essential and "UP2" not in essential
    assert "R1" in essential


def test_an_impossible_task_has_no_essential_reactions(task_model: Model) -> None:
    """Nothing is essential to something that cannot be done at all."""
    impossible = MetabolicTask(
        id="impossible",
        outputs={"atp_c": (1.0, 1000.0)},
    )

    assert essential_reactions_for_task(task_model, impossible) == set()


def test_essential_reactions_union_over_tasks(task_model: Model) -> None:
    """The union is what becomes the core set, and should-fail tasks are out."""
    tasks = [
        _atp_task(),
        MetabolicTask(
            id="waste",
            inputs={"s_e": (0.0, 1000.0)},
            outputs={"w_c": (1.0, 1000.0)},
        ),
        MetabolicTask(id="control", equations=(ATP_HYDROLYSIS,), should_fail=True),
    ]

    union, per_task = essential_reactions_for_tasks(task_model, tasks)

    assert "control" not in per_task, "should-fail tasks must not force anything in"
    assert per_task["atp"] == {"UP", "R1"}
    assert union == set().union(*per_task.values())
    assert "UP" in union


def test_check_tasks_reports_each_in_order(task_model: Model) -> None:
    """The report lines up with the tasks given."""
    tasks = [
        _atp_task(),
        MetabolicTask("control", equations=(ATP_HYDROLYSIS,), should_fail=True),
    ]

    results = check_tasks(task_model, tasks)

    assert [r.task.id for r in results] == ["atp", "control"]
    assert all(r.ok for r in results)
