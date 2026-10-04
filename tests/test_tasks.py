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


def test_free_inputs_and_outputs_loosen_every_task(task_model: Model) -> None:
    """Published lists omit what they take for granted, so it can be supplied.

    A task is written against a convention: common waste can leave, oxygen
    is available. Where the convention is not written down, the task is
    infeasible as read, and these two options are how the caller restores
    it.
    """
    import tempfile
    from pathlib import Path

    from cobra import Metabolite
    from cobra import Model as CobraModel
    from cobra import Reaction

    from spectra_cobra.tasks import parse_task_list

    # The parser resolves by metabolite *name* and compartment, so the model
    # it binds to has to carry both.
    named = CobraModel("named")
    substrate = Metabolite("s", name="substrate", compartment="e")
    waste = Metabolite("w", name="waste", compartment="c")
    named.add_metabolites([substrate, waste])
    convert = Reaction("CONV", lower_bound=0.0, upper_bound=1000.0)
    named.add_reactions([convert])
    convert.add_metabolites({substrate: -1.0, waste: 1.0})
    named.solver = task_model.solver.interface.__name__.rsplit(".", 1)[-1]

    header = (
        "\tID\tDESCRIPTION\tSHOULD FAIL\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB"
        "\tEQU\tEQU LB\tEQU UB\n"
    )
    # Consume the substrate; the waste it makes is deliberately not listed.
    row = "\tT\tconsume the substrate\t\tsubstrate[e]\t1\t1\t\t\t\t\t\t\n"
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "tasks.txt"
        path.write_text(header + row)

        bare, _ = parse_task_list(str(path), named)
        loosened, _ = parse_task_list(str(path), named, free_outputs=("waste[c]",))

    assert bare and loosened
    # Consuming the substrate makes waste, and the row gives it nowhere to go.
    assert not check_task(named, bare[0]).feasible
    assert check_task(named, loosened[0]).feasible


def test_an_equation_keeps_the_proton_together(task_model: Model) -> None:
    """Terms split on " + ", so "H+[c]" must survive as one reference.

    Splitting on a bare "+" shears the proton into "H" and "[c]", which cost
    26 of the 69 rows of the Human-GEM list before it was fixed.
    """
    from cobra import Metabolite
    from cobra import Model as CobraModel

    from spectra_cobra.tasks import _metabolite_index, _parse_equation

    model = CobraModel("protons")
    model.add_metabolites(
        [
            Metabolite("atp", name="ATP", compartment="c"),
            Metabolite("h2o", name="H2O", compartment="c"),
            Metabolite("adp", name="ADP", compartment="c"),
            Metabolite("h", name="H+", compartment="c"),
        ]
    )
    index = _metabolite_index(model)

    stoichiometry, unresolved = _parse_equation(
        "ATP[c] + H2O[c] => ADP[c] + H+[c]", index
    )

    assert unresolved == [], f"the proton was sheared apart: {unresolved}"
    assert stoichiometry == {"atp": -1.0, "h2o": -1.0, "adp": 1.0, "h": 1.0}


def test_a_metabolite_in_both_lists_may_flow_either_way(solver: str) -> None:
    """Naming a metabolite twice says it may go both ways, and it does.

    The ambiguity belongs to the published file format, which has one bound
    column per direction and no way to say which was meant; the parser
    resolves it there. A task built directly says what it says, so with an
    input of (1, 1) and an output of (0, 1000) the forced unit of supply
    takes the top of the range from 1000 down to 999: the drain has to
    carry that unit as well as everything the network made.
    """
    from cobra import Metabolite
    from cobra import Model as CobraModel
    from cobra import Reaction

    from spectra_cobra.tasks import task_constraints

    def build() -> CobraModel:
        model = CobraModel("both")
        made = Metabolite("M", compartment="c")
        fed = Metabolite("S", compartment="c")
        model.add_metabolites([made, fed])
        convert = Reaction("R", lower_bound=-1000.0, upper_bound=1000.0)
        source = Reaction("SRC", lower_bound=-1000.0, upper_bound=1000.0)
        model.add_reactions([convert, source])
        convert.add_metabolites({fed: -1.0, made: 1.0})
        source.add_metabolites({fed: 1.0})
        model.solver = solver
        return model

    def net_range(task) -> tuple:
        model = build()
        with task_constraints(model, task) as constrained:
            built = {
                r.id for r in constrained.reactions if r.id.startswith("spectra_task_")
            }
            span = []
            for direction in ("max", "min"):
                constrained.objective = constrained.reactions.R
                constrained.objective.direction = direction
                span.append(constrained.slim_optimize())
        return built, (round(span[1], 6), round(span[0], 6))

    feed = {"S": (0.0, 1000.0)}
    output_only, span_out = net_range(
        MetabolicTask("t", inputs=feed, outputs={"M": (0.0, 1000.0)})
    )
    both, span_both = net_range(
        MetabolicTask(
            "t", inputs={**feed, "M": (1.0, 1.0)}, outputs={"M": (0.0, 1000.0)}
        )
    )

    assert span_out == (0.0, 1000.0)
    # The lower end stays at 0 because the supply of S cannot run backwards,
    # not because the input declaration was dropped; the upper end is what
    # shows it was honoured.
    assert span_both == (0.0, 999.0), "both declarations should be honoured"
    assert len(both) == len(output_only) + 1


def test_the_parser_resolves_a_row_that_says_both(task_model: Model) -> None:
    """A metabolite in both columns follows the published format's rule.

    The output's upper bound always wins and the input's lower bound is
    discarded, but the output's lower bound wins only when it is positive.
    So an output that merely *permits* production leaves the input's
    allowance to consume intact and the metabolite may still flow either
    way; an output that *requires* production replaces the input outright.
    Requiring both at once is contradictory and raises.
    """
    import tempfile
    from pathlib import Path

    from cobra import Metabolite
    from cobra import Model as CobraModel
    from cobra import Reaction

    from spectra_cobra.tasks import parse_task_list

    model = CobraModel("named")
    made = Metabolite("m", name="made", compartment="c")
    fed = Metabolite("f", name="fed", compartment="c")
    model.add_metabolites([made, fed])
    convert = Reaction("R", lower_bound=0.0, upper_bound=1000.0)
    model.add_reactions([convert])
    convert.add_metabolites({fed: -1.0, made: 1.0})

    header = (
        "\tID\tDESCRIPTION\tSHOULD FAIL\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB"
        "\tEQU\tEQU LB\tEQU UB\n"
    )
    # OUT LB = 0, so the output only permits production: the input's
    # allowance to consume up to 1 survives, its requirement does not.
    permits = "\tT\tboth ways\t\tfed[c];made[c]\t1\t1\tmade[c]\t0\t1000\t\t\t\n"
    # OUT LB = 5, so the output requires production and replaces the input.
    requires = "\tT\tboth ways\t\tfed[c];made[c]\t0\t1\tmade[c]\t5\t1000\t\t\t\n"
    # Both lower bounds positive is contradictory.
    clashes = "\tT\tboth ways\t\tfed[c];made[c]\t1\t1\tmade[c]\t5\t1000\t\t\t\n"

    def parsed(row, **kwargs):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "tasks.txt"
            path.write_text(header + row)
            return parse_task_list(str(path), model, **kwargs)[0]

    permitting = parsed(permits)[0]
    assert permitting.inputs["m"] == (0.0, 1.0), "the allowance should survive"
    assert permitting.outputs["m"] == (0.0, 1000.0)

    requiring = parsed(requires)[0]
    assert "m" not in requiring.inputs, "a required output replaces the input"

    with pytest.raises(SpectraError, match="cannot both hold"):
        parsed(clashes)

    # free_inputs is the caller's word, not the file's, so it is untouched.
    assert "m" in parsed(requires, free_inputs=("made[c]",))[0].inputs
