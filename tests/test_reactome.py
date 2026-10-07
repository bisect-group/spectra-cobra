"""Test cutting a model down to the smallest network that still works."""

import logging

import pytest
from cobra import Metabolite, Model, Reaction

from spectra_cobra import (
    MetabolicTask,
    SpectraError,
    TaskEquation,
    check_tasks,
    minimal_reactome,
)


def _add(model: Model, rxn_id: str, stoichiometry: dict, bounds=(0.0, 1000.0)) -> None:
    """Add one reaction to a model.

    Parameters
    ----------
    model : cobra.Model
        The model to add to.
    rxn_id : str
        The reaction's identifier.
    stoichiometry : dict of {str: float}
        The coefficients, keyed by metabolite identifier.
    bounds : tuple of (float, float), optional
        The flux bounds (default ``(0, 1000)``).

    """
    reaction = Reaction(rxn_id, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([reaction])
    reaction.add_metabolites(
        {
            model.metabolites.get_by_id(met): coeff
            for met, coeff in stoichiometry.items()
        }
    )


@pytest.fixture(scope="function")
def cell(solver: str) -> Model:
    """Return a cell that can live on either of two substrates.

    Glucose reaches the growth precursor ``x_c`` either in one step
    (``T_GLC``) or in two (``T_GLC_LONG_1`` and ``T_GLC_LONG_2``); acetate
    reaches it in one (``T_AC``); and ``P_AC`` runs the other way, turning
    the precursor back into acetate so that secretion can be required of
    it. So the smallest network that grows is three reactions, and *which*
    three is decided entirely by the medium.

    Returns
    -------
    cobra.Model
        The toy cell.

    """
    model = Model("cell")
    model.add_metabolites(
        [
            Metabolite("glc_e", compartment="e"),
            Metabolite("ac_e", compartment="e"),
            Metabolite("x_c", compartment="c"),
            Metabolite("y_c", compartment="c"),
        ]
    )
    _add(model, "T_GLC", {"glc_e": -1.0, "x_c": 1.0})
    _add(model, "T_GLC_LONG_1", {"glc_e": -1.0, "y_c": 1.0})
    _add(model, "T_GLC_LONG_2", {"y_c": -1.0, "x_c": 1.0})
    _add(model, "T_AC", {"ac_e": -1.0, "x_c": 1.0})
    _add(model, "P_AC", {"x_c": -1.0, "ac_e": 1.0})
    _add(model, "BIO", {"x_c": -1.0})
    model.add_boundary(model.metabolites.get_by_id("glc_e"), type="exchange")
    model.add_boundary(model.metabolites.get_by_id("ac_e"), type="exchange")
    model.objective = model.reactions.BIO
    model.solver = solver
    model.tolerance = 1e-9
    return model


GLUCOSE = {"EX_glc_e": -10.0}
ACETATE = {"EX_ac_e": -10.0}


def test_the_medium_decides_which_reactions_survive(cell) -> None:
    """Burgard's finding in miniature: minimal for an environment."""
    on_glucose = minimal_reactome(cell, medium=GLUCOSE)
    on_acetate = minimal_reactome(cell, medium=ACETATE)

    assert set(on_glucose.kept) == {"EX_glc_e", "T_GLC", "BIO"}
    assert set(on_acetate.kept) == {"EX_ac_e", "T_AC", "BIO"}


def test_the_count_is_minimised_not_the_flux(cell) -> None:
    """The one-step route wins over the two-step one on reaction count."""
    result = minimal_reactome(cell, medium=GLUCOSE)

    assert "T_GLC" in result.kept
    assert "T_GLC_LONG_1" not in result.kept
    assert "T_GLC_LONG_2" not in result.kept


def test_a_weight_can_make_the_longer_route_the_cheaper_one(cell) -> None:
    """Weights are what a reaction costs to keep, not how much flux it gets."""
    result = minimal_reactome(cell, medium=GLUCOSE, weights={"T_GLC": 5.0})

    assert "T_GLC" not in result.kept
    assert {"T_GLC_LONG_1", "T_GLC_LONG_2"} <= set(result.kept)


def test_growth_is_preserved_in_full_by_default(cell) -> None:
    """The default demands every bit of the growth the medium allows."""
    result = minimal_reactome(cell, medium=GLUCOSE)

    with cell:
        for reaction in cell.boundary:
            reaction.lower_bound = 0.0
        cell.reactions.EX_glc_e.lower_bound = -10.0
        full = cell.slim_optimize()
    assert result.required_growth["growth"] == pytest.approx(full)
    assert result.growth["growth"] == pytest.approx(full)
    assert result.unsatisfied == ()


def test_a_lower_requirement_does_not_need_a_smaller_network_here(cell) -> None:
    """Three reactions is the floor however little growth is asked for."""
    result = minimal_reactome(cell, medium=GLUCOSE, growth_fraction=0.1)

    assert len(result.kept) == 3
    assert result.required_growth["growth"] == pytest.approx(1.0)


def test_an_absolute_growth_rate_can_be_demanded_instead(cell) -> None:
    """min_growth is a bound on the biomass reaction, not a fraction."""
    result = minimal_reactome(cell, medium=GLUCOSE, min_growth=4.0)

    assert result.required_growth["growth"] == pytest.approx(4.0)
    assert result.growth["growth"] >= 4.0


def test_a_required_product_keeps_the_reactions_that_make_it(cell) -> None:
    """Secretion is a requirement like any other, and costs reactions."""
    without = minimal_reactome(cell, medium=GLUCOSE, growth_fraction=0.5)
    with_acetate = minimal_reactome(
        cell, medium=GLUCOSE, growth_fraction=0.5, products={"EX_ac_e": 1.0}
    )

    assert "P_AC" not in without.kept
    assert {"P_AC", "EX_ac_e"} <= set(with_acetate.kept)
    assert with_acetate.production["growth"]["EX_ac_e"] >= 1.0
    assert with_acetate.unsatisfied == ()


def test_a_product_named_as_a_metabolite_uses_the_exchange_it_has(cell) -> None:
    """Naming the metabolite and naming its exchange are the same request."""
    by_metabolite = minimal_reactome(
        cell, medium=GLUCOSE, growth_fraction=0.5, products={"ac_e": 1.0}
    )
    by_reaction = minimal_reactome(
        cell, medium=GLUCOSE, growth_fraction=0.5, products={"EX_ac_e": 1.0}
    )

    assert set(by_metabolite.kept) == set(by_reaction.kept)
    assert list(by_metabolite.production["growth"]) == ["EX_ac_e"]
    assert "EX_ac_e" in by_metabolite.kept


def test_nothing_is_added_to_the_network_to_measure_a_product(cell) -> None:
    """The result has to stay a subnetwork of the model handed in."""
    before = {rxn.id for rxn in cell.reactions}

    result = minimal_reactome(
        cell, medium=GLUCOSE, growth_fraction=0.5, products={"ac_e": 1.0}
    )

    assert {rxn.id for rxn in cell.reactions} == before
    assert set(result.kept) <= before
    assert not any(rxn.id.startswith("DM_") for rxn in result.model.reactions)


def test_an_unexchanged_metabolite_is_refused_rather_than_drained(cell) -> None:
    """A drain added here would change the network being minimised.

    ``y_c`` sits between the two halves of the long glucose route and no
    boundary reaction touches it. Adding one would make the answer the
    minimal subnetwork of a model the caller does not have, feasible only
    because of the reaction added to measure it.
    """
    with pytest.raises(SpectraError, match="no exchange reaction"):
        minimal_reactome(
            cell, medium=GLUCOSE, growth_fraction=0.5, products={"y_c": 1.0}
        )


def test_a_metabolite_with_two_boundary_reactions_is_ambiguous(cell) -> None:
    """Which one the production is measured through would be a coin toss."""
    cell.add_boundary(cell.metabolites.get_by_id("ac_e"), type="demand")

    with pytest.raises(SpectraError, match="ambiguous"):
        minimal_reactome(
            cell, medium=GLUCOSE, growth_fraction=0.5, products={"ac_e": 1.0}
        )


def test_an_exchange_written_as_a_source_is_refused(cell) -> None:
    """Flux through it is uptake, so a lower bound would require the opposite."""
    cell.add_metabolites([Metabolite("z_e", compartment="e")])
    _add(cell, "SRC_z", {"z_e": 1.0})

    with pytest.raises(SpectraError, match="written as a source"):
        minimal_reactome(
            cell, medium=GLUCOSE, growth_fraction=0.5, products={"z_e": 1.0}
        )


def test_growth_and_production_have_to_hold_at_once(cell) -> None:
    """Both are the same condition, so a network that can only alternate fails."""
    # All ten units of glucose are needed for full growth, so none is left
    # to become acetate -- though either on its own is perfectly reachable.
    with pytest.raises(SpectraError, match="cannot meet condition"):
        minimal_reactome(
            cell, medium=GLUCOSE, growth_fraction=1.0, products={"EX_ac_e": 1.0}
        )


def test_an_impossible_pair_says_which_half_was_reachable(cell) -> None:
    """An infeasible MILP names nothing; this has to name the requirements."""
    with pytest.raises(SpectraError) as raised:
        minimal_reactome(
            cell, medium=GLUCOSE, growth_fraction=1.0, products={"EX_ac_e": 1.0}
        )

    message = str(raised.value)
    assert "BIO at 10" in message
    assert "EX_ac_e at 1" in message
    assert message.count("(reachable)") == 2, "each alone, but not together"


def test_a_kept_reaction_has_to_carry_flux_not_merely_be_present(cell) -> None:
    """Forcing T_AC drags in the only thing that can feed it.

    On glucose nothing takes up acetate, so the one source of ``ac_e`` is
    ``P_AC`` running the other way. Retaining ``T_AC`` without requiring
    flux would leave it dead and ``P_AC`` out; requiring flux cannot.
    """
    result = minimal_reactome(cell, medium=GLUCOSE, keep_reactions=["T_AC"])

    assert set(result.kept) == {"EX_glc_e", "T_GLC", "BIO", "T_AC", "P_AC"}
    assert result.growth["growth"] == pytest.approx(10.0), "and it still has to grow"


def test_forcing_half_a_route_brings_in_the_other_half(cell) -> None:
    """And drops the one-step route it was competing with."""
    result = minimal_reactome(cell, medium=GLUCOSE, keep_reactions=["T_GLC_LONG_1"])

    assert set(result.kept) == {
        "EX_glc_e",
        "T_GLC_LONG_1",
        "T_GLC_LONG_2",
        "BIO",
    }
    assert "T_GLC" not in result.kept


def test_a_kept_reaction_that_cannot_carry_flux_is_reported(cell) -> None:
    """Retained anyway, but you are told nothing is making it work."""
    cell.add_metabolites([Metabolite("dead_c", compartment="c")])
    _add(cell, "DEAD", {"x_c": -1.0, "dead_c": 1.0})

    result = minimal_reactome(cell, medium=GLUCOSE, keep_reactions=["DEAD"])

    assert result.kept_but_blocked == {"growth": ("DEAD",)}
    assert "DEAD" in result.kept, "still retained -- it was asked for"
    assert "unable to carry flux" in result.summary()


def test_kept_reactions_are_not_forced_inside_a_task(cell) -> None:
    """A task closes every exchange, so forcing one there is impossible.

    This is why the flux requirement applies to the growth conditions
    only. Under a rule that forced ``keep_reactions`` in every condition,
    keeping any boundary reaction alongside any task would be infeasible
    by construction.
    """
    on_acetate = MetabolicTask(
        id="precursor_from_acetate",
        inputs={"ac_e": (0.0, 10.0)},
        equations=[_drain("x_c")],
    )

    result = minimal_reactome(
        cell, medium=GLUCOSE, keep_reactions=["EX_glc_e"], tasks=[on_acetate]
    )

    assert "EX_glc_e" in result.kept
    assert result.unsatisfied == ()
    assert all(outcome.ok for outcome in check_tasks(result.model, [on_acetate]))


def test_the_model_handed_in_is_not_modified(cell) -> None:
    """The result is a new model; nothing is removed from the caller's."""
    before = {rxn.id: rxn.bounds for rxn in cell.reactions}

    result = minimal_reactome(cell, medium=GLUCOSE)

    assert {rxn.id: rxn.bounds for rxn in cell.reactions} == before
    assert len(result.model.reactions) == len(result.kept)


def test_alternative_minimal_reactomes_are_distinct(cell) -> None:
    """Each solution is forbidden from recurring before the next solve."""
    results = minimal_reactome(cell, medium=GLUCOSE, growth_fraction=0.1, n_solutions=2)

    assert len(results) == 2
    assert {frozenset(r.kept) for r in results} == {
        frozenset({"EX_glc_e", "T_GLC", "BIO"}),
        frozenset({"EX_glc_e", "T_GLC_LONG_1", "T_GLC_LONG_2", "BIO"}),
    }
    assert all(r.growth["growth"] >= r.required_growth["growth"] for r in results)


def test_asking_for_more_alternatives_than_exist_returns_what_there_is(
    cell, caplog
) -> None:
    """Two routes to the precursor, so there is no third answer to give."""
    with caplog.at_level(logging.INFO, logger="spectra_cobra.extraction"):
        results = minimal_reactome(
            cell, medium=GLUCOSE, growth_fraction=0.1, n_solutions=5
        )

    assert len(results) == 2
    assert "No further alternative solution exists" in caplog.text


def test_the_summary_says_what_happened(cell) -> None:
    """A result is useless if you cannot read it at a glance."""
    result = minimal_reactome(cell, medium=GLUCOSE)

    assert "3 of 8 reactions" in result.summary()
    assert "growth" in result.summary()


# --------------------------------------------------------------------------
# Several media, which are separate conditions just as tasks are.
# --------------------------------------------------------------------------


def test_several_media_must_all_work(cell) -> None:
    """The network has to grow in every one of them, not in the easiest."""
    result = minimal_reactome(cell, medium={"glucose": GLUCOSE, "acetate": ACETATE})

    assert set(result.kept) == {"EX_glc_e", "T_GLC", "EX_ac_e", "T_AC", "BIO"}
    assert result.conditions == ("glucose", "acetate")
    assert set(result.growth) == {"glucose", "acetate"}
    assert result.unsatisfied == ()


def test_each_medium_gets_its_own_growth_requirement(cell) -> None:
    """A fraction is of what that medium allows, not of some global rate."""
    halved = dict(GLUCOSE, EX_glc_e=-4.0)

    result = minimal_reactome(
        cell, medium={"rich": GLUCOSE, "poor": halved}, growth_fraction=1.0
    )

    assert result.required_growth["rich"] == pytest.approx(10.0)
    assert result.required_growth["poor"] == pytest.approx(4.0)
    assert result.growth["poor"] == pytest.approx(4.0)


def test_media_can_be_given_without_labels(cell) -> None:
    """A plain sequence is labelled for you."""
    result = minimal_reactome(cell, medium=[GLUCOSE, ACETATE])

    assert result.conditions == ("medium_0", "medium_1")
    assert len(result.kept) == 5


def test_one_medium_in_a_sequence_is_still_one_medium(cell) -> None:
    """No spurious label for a list of one."""
    result = minimal_reactome(cell, medium=[GLUCOSE])

    assert result.conditions == ("growth",)
    assert set(result.kept) == {"EX_glc_e", "T_GLC", "BIO"}


def test_a_medium_mixing_bounds_and_media_is_refused(cell) -> None:
    """Whether it is one medium or several would otherwise be a guess."""
    with pytest.raises(SpectraError, match="mixes"):
        minimal_reactome(cell, medium={"EX_glc_e": -10.0, "other": {"EX_ac_e": -10.0}})


# --------------------------------------------------------------------------
# Preprocessing, which must change the cost and not the answer.
# --------------------------------------------------------------------------


def test_preprocessing_does_not_change_the_answer(cell) -> None:
    """Both passes are exact, so only the work should change.

    What is guaranteed is the optimum, not the last reaction of the
    readout: the trace-flux rule behind ``inclusion_cutoff`` can sweep in
    a reaction whose binary is off, and a differently shaped relaxation
    leaves traces in different places. On a model this small there is
    nowhere for a trace to hide, so the sets match exactly.
    """
    with_it = minimal_reactome(cell, medium=GLUCOSE, preprocess=True)
    without = minimal_reactome(cell, medium=GLUCOSE, preprocess=False)

    assert set(with_it.kept) == set(without.kept)
    assert with_it.essential, "something has to be essential here"
    assert without.essential == ()


def test_essential_means_the_condition_fails_without_it(cell) -> None:
    """Not merely 'in the answer' -- there must be no alternative at all."""
    result = minimal_reactome(cell, medium=GLUCOSE)

    assert set(result.essential) <= set(result.kept)
    assert {"BIO", "EX_glc_e"} <= set(result.essential)
    assert "T_GLC" not in result.essential, "the two-step route is an alternative"


def test_a_blocked_reaction_is_deleted_and_still_reported(cell) -> None:
    """It was never a candidate, but it is still not in the answer."""
    cell.add_metabolites([Metabolite("dead_c", compartment="c")])
    _add(cell, "DEAD", {"x_c": -1.0, "dead_c": 1.0})

    result = minimal_reactome(cell, medium=GLUCOSE)

    assert "DEAD" not in result.model.reactions
    assert "DEAD" in result.removed
    assert len(result.kept) + len(result.removed) == len(cell.reactions)


def test_kept_reactions_survive_preprocessing(cell) -> None:
    """A reaction the caller pinned is never second-guessed."""
    result = minimal_reactome(cell, medium=GLUCOSE, keep_reactions=["P_AC"])

    assert "P_AC" in result.kept
    assert "P_AC" not in result.essential, "it is kept by request, not by need"


def _drain(metabolite: str, rate: float = 1.0) -> TaskEquation:
    """Return an equation that forces a metabolite to be made.

    Parameters
    ----------
    metabolite : str
        What has to be produced.
    rate : float, optional
        How fast (default 1).

    Returns
    -------
    TaskEquation
        The forcing equation.

    """
    return TaskEquation({metabolite: -1.0}, (rate, 1000.0))


def test_a_task_keeps_reactions_growth_alone_would_drop(cell) -> None:
    """The task runs on acetate, so the acetate transporter must survive."""
    on_acetate = MetabolicTask(
        id="precursor_from_acetate",
        inputs={"ac_e": (0.0, 10.0)},
        equations=[_drain("x_c")],
    )

    result = minimal_reactome(cell, medium=GLUCOSE, tasks=[on_acetate])

    assert "T_AC" in result.kept, "nothing about growing on glucose needs this"
    assert set(result.kept) == {"EX_glc_e", "T_GLC", "BIO", "T_AC"}
    assert result.conditions == ("growth", "precursor_from_acetate")
    assert result.unsatisfied == ()


def test_the_result_really_performs_the_tasks(cell) -> None:
    """Verified against the task checker, not just against the solver."""
    on_acetate = MetabolicTask(
        id="precursor_from_acetate",
        inputs={"ac_e": (0.0, 10.0)},
        equations=[_drain("x_c")],
    )

    result = minimal_reactome(cell, medium=GLUCOSE, tasks=[on_acetate])

    assert all(outcome.ok for outcome in check_tasks(result.model, [on_acetate]))


def test_a_task_the_model_cannot_do_is_an_error_not_an_infeasibility(cell) -> None:
    """Removing reactions cannot make a task work, so say so plainly."""
    impossible = MetabolicTask(
        id="precursor_from_nothing",
        inputs={},
        equations=[_drain("x_c")],
    )

    with pytest.raises(SpectraError, match="cannot perform"):
        minimal_reactome(cell, medium=GLUCOSE, tasks=[impossible])


@pytest.fixture(scope="function")
def two_routes(solver: str) -> Model:
    """Return a network where sharing beats two private routes.

    Each substrate reaches ``out_c`` either by its own three-reaction
    private route, or by a one-reaction entry into a three-reaction hub
    the other substrate can use too. On its own a substrate prefers its
    private route, three reactions against four. Both together prefer the
    hub, five reactions against six. So solving the two conditions apart
    and taking the union gives six, and solving them jointly gives five.

    Returns
    -------
    cobra.Model
        The toy network.

    """
    model = Model("two_routes")
    model.add_metabolites(
        [Metabolite(f"{name}", compartment="c") for name in ("m1_c", "m2_c", "m3_c")]
        + [Metabolite(f"{name}", compartment="c") for name in ("q1_c", "q2_c")]
        + [Metabolite(f"{name}", compartment="c") for name in ("r1_c", "r2_c")]
        + [Metabolite("out_c", compartment="c")]
        + [Metabolite("s1_e", compartment="e"), Metabolite("s2_e", compartment="e")]
    )
    _add(model, "E1", {"s1_e": -1.0, "m1_c": 1.0})
    _add(model, "E2", {"s2_e": -1.0, "m1_c": 1.0})
    _add(model, "H1", {"m1_c": -1.0, "m2_c": 1.0})
    _add(model, "H2", {"m2_c": -1.0, "m3_c": 1.0})
    _add(model, "H3", {"m3_c": -1.0, "out_c": 1.0})
    _add(model, "Q1a", {"s1_e": -1.0, "q1_c": 1.0})
    _add(model, "Q1b", {"q1_c": -1.0, "q2_c": 1.0})
    _add(model, "Q1c", {"q2_c": -1.0, "out_c": 1.0})
    _add(model, "Q2a", {"s2_e": -1.0, "r1_c": 1.0})
    _add(model, "Q2b", {"r1_c": -1.0, "r2_c": 1.0})
    _add(model, "Q2c", {"r2_c": -1.0, "out_c": 1.0})
    model.add_boundary(model.metabolites.get_by_id("s1_e"), type="exchange")
    model.add_boundary(model.metabolites.get_by_id("s2_e"), type="exchange")
    model.solver = solver
    model.tolerance = 1e-9
    return model


def _from(substrate: str, label: str) -> MetabolicTask:
    """Return a task making ``out_c`` from one substrate.

    Parameters
    ----------
    substrate : str
        The extracellular metabolite to feed in.
    label : str
        The task's identifier.

    Returns
    -------
    MetabolicTask
        The task.

    """
    return MetabolicTask(
        id=label,
        inputs={substrate: (0.0, 10.0)},
        equations=[_drain("out_c")],
    )


def test_one_condition_at_a_time_takes_the_private_route(two_routes) -> None:
    """Three reactions beats four, when there is nobody to share with."""
    alone = minimal_reactome(
        two_routes, growth_fraction=0.0, tasks=[_from("s1_e", "a")]
    )

    assert set(alone.kept) == {"Q1a", "Q1b", "Q1c"}


def test_solving_the_conditions_jointly_beats_their_union(two_routes) -> None:
    """The whole reason the conditions share one set of binaries."""
    first = _from("s1_e", "a")
    second = _from("s2_e", "b")

    joint = minimal_reactome(two_routes, growth_fraction=0.0, tasks=[first, second])
    union = set(
        minimal_reactome(two_routes, growth_fraction=0.0, tasks=[first]).kept
    ) | set(minimal_reactome(two_routes, growth_fraction=0.0, tasks=[second]).kept)

    assert set(joint.kept) == {"E1", "E2", "H1", "H2", "H3"}
    assert len(joint.kept) == 5
    assert len(union) == 6
    assert all(outcome.ok for outcome in check_tasks(joint.model, [first, second]))


# --------------------------------------------------------------------------
# Arguments that cannot mean anything.
# --------------------------------------------------------------------------


def test_two_ways_of_asking_for_growth_at_once_is_refused(cell) -> None:
    """One is a rate and the other a fraction; both leaves it ambiguous."""
    with pytest.raises(SpectraError, match="not both"):
        minimal_reactome(cell, medium=GLUCOSE, min_growth=1.0, growth_fraction=0.5)


def test_a_medium_naming_an_absent_reaction_is_refused(cell) -> None:
    """Silently ignoring it would silently change the environment."""
    with pytest.raises(SpectraError, match="does not"):
        minimal_reactome(cell, medium={"EX_nope_e": -10.0})


def test_a_product_that_is_neither_reaction_nor_metabolite_is_refused(cell) -> None:
    """There would be nothing to put a bound on."""
    with pytest.raises(SpectraError, match="neither"):
        minimal_reactome(cell, medium=GLUCOSE, products={"nope": 1.0})


def test_a_production_rate_of_zero_is_refused(cell) -> None:
    """Zero flux satisfies it, so it is not a requirement at all."""
    with pytest.raises(SpectraError, match="requires nothing"):
        minimal_reactome(cell, medium=GLUCOSE, products={"EX_ac_e": 0.0})


def test_an_unknown_reaction_to_keep_is_refused(cell) -> None:
    """A typo here would quietly drop the protection it was asking for."""
    with pytest.raises(SpectraError, match="keep_reactions"):
        minimal_reactome(cell, medium=GLUCOSE, keep_reactions=["NOPE"])


def test_a_negative_weight_is_refused(cell) -> None:
    """A negative cost would pay the solver to keep reactions."""
    with pytest.raises(SpectraError, match="non-negative"):
        minimal_reactome(cell, medium=GLUCOSE, weights={"T_AC": -1.0})


def test_a_fraction_of_growth_that_cannot_happen_is_refused(cell) -> None:
    """There is no rate to take a fraction of, so say that rather than fail."""
    with pytest.raises(SpectraError, match="no rate to take a fraction of"):
        minimal_reactome(cell, medium={}, growth_fraction=0.5)


def test_asking_for_nothing_is_refused(cell) -> None:
    """The smallest network meeting no requirement is the empty one."""
    with pytest.raises(SpectraError, match="Nothing was asked"):
        minimal_reactome(cell, growth_fraction=0.0)


def test_an_ambiguous_objective_must_be_resolved_when_growth_matters(cell) -> None:
    """Which reaction is growth cannot be guessed from two of them."""
    cell.objective = {cell.reactions.BIO: 1.0, cell.reactions.P_AC: 1.0}

    with pytest.raises(SpectraError, match="ambiguous"):
        minimal_reactome(cell, medium=GLUCOSE, growth_fraction=0.5)
