"""Test model extraction with each network inference formulation."""

from typing import List

import pytest
from cobra import Model

from spectra_cobra import (
    GROWTH_OPTIM,
    MIN_NET_LP,
    MIN_NET_MILP,
    PATHWAY_EXCLUSION,
    TRADE_OFF,
    SpectraError,
    SpectraInfeasibleCoreError,
    spectra_ccme,
    spectra_me,
)

#: Every formulation that takes non-negative weights, i.e. all but tradeOff.
COST_PROBLEM_TYPES = (MIN_NET_LP, MIN_NET_MILP)


def _ids(model: Model) -> set:
    """Return a model's reaction identifiers.

    Parameters
    ----------
    model : cobra.Model
        The model to inspect.

    Returns
    -------
    set of str
        The reaction identifiers.

    """
    return {rxn.id for rxn in model.reactions}


@pytest.mark.parametrize("problem_type", COST_PROBLEM_TYPES)
def test_extracts_the_pathway_through_the_core(
    toy_model: Model, problem_type: str
) -> None:
    """A core reaction pulls in exactly the pathway that feeds it.

    ``R3`` sits on the A-B-C-D pathway, so keeping it requires R1, R2, R4 and
    R5 and nothing from the other two pathways.
    """
    extracted = spectra_me(toy_model, ["R3"], problem_type=problem_type, seed=0)
    assert _ids(extracted) == {"R1", "R2", "R3", "R4", "R5"}


def test_core_reaction_is_always_kept(toy_model: Model) -> None:
    """Each pathway's core reaction appears in its own extracted model."""
    for core, expected in (
        ("R3", {"R1", "R2", "R3", "R4", "R5"}),
        ("R7", {"R6", "R7", "R8", "R5"}),
        ("R10", {"R9", "R10", "R11", "R5"}),
    ):
        extracted = spectra_me(toy_model, [core], seed=0)
        assert core in _ids(extracted)
        assert _ids(extracted) == expected


def test_weights_steer_which_pathway_is_kept(toy_model: Model) -> None:
    """Making one pathway expensive pushes the solver onto another.

    With ``R5`` as the only core reaction all three pathways can satisfy it,
    so the weights decide. Penalising a pathway's reactions heavily should
    keep it out of the extracted model.
    """
    cheap, expensive = 1.0, 1000.0
    weights = {r.id: cheap for r in toy_model.reactions}
    for rxn_id in ("R1", "R2", "R3", "R4", "R6", "R7", "R8"):
        weights[rxn_id] = expensive

    extracted = spectra_me(
        toy_model, ["R5"], weights=weights, problem_type=MIN_NET_MILP, seed=0
    )
    # The G-H pathway is the only cheap one left.
    assert {"R9", "R10", "R11"} <= _ids(extracted)
    assert not {"R2", "R3", "R7"} & _ids(extracted)


def test_trade_off_accepts_negative_weights(toy_model: Model) -> None:
    """A negatively weighted reaction is pushed out of the model.

    ``tradeOff`` maximises the weighted count of kept reactions, so a
    reaction with a negative weight is only kept if the network cannot do
    without it.
    """
    weights = {r.id: 1.0 for r in toy_model.reactions}
    for rxn_id in ("R7", "R8", "R10", "R11"):
        weights[rxn_id] = -5.0

    extracted = spectra_me(
        toy_model, ["R3"], weights=weights, problem_type=TRADE_OFF, seed=0
    )
    assert {"R1", "R2", "R3", "R4", "R5"} <= _ids(extracted)
    assert not {"R7", "R8", "R10", "R11"} & _ids(extracted)


def test_growth_optim_keeps_a_growing_network(toy_model: Model) -> None:
    """growthOptim retains a route to the objective reaction."""
    extracted = spectra_me(toy_model, ["R3"], problem_type=GROWTH_OPTIM, seed=0)
    assert "R5" in _ids(extracted)
    assert {"R1", "R2", "R3", "R4"} <= _ids(extracted)


def test_growth_optim_needs_an_objective(toy_model: Model) -> None:
    """growthOptim without an objective reaction is an error."""
    toy_model.objective = toy_model.problem.Objective(0)
    with pytest.raises(SpectraError, match="objective"):
        spectra_me(toy_model, ["R3"], problem_type=GROWTH_OPTIM, seed=0)


def test_blocked_core_reaction_raises(blocked_model: Model) -> None:
    """spectra_me refuses a core reaction that cannot carry flux."""
    with pytest.raises(SpectraInfeasibleCoreError, match="spectra_ccme"):
        spectra_me(blocked_model, ["R3"], seed=0)


def test_ccme_reports_and_drops_blocked_core(blocked_model: Model) -> None:
    """spectra_ccme drops a blocked core reaction and names it."""
    extracted, blocked_core = spectra_ccme(blocked_model, ["R1", "R3"], seed=0)
    assert blocked_core == ["R3"]
    assert "R3" not in _ids(extracted)
    assert "R1" in _ids(extracted)


def test_ccme_matches_me_on_a_consistent_model(toy_model: Model) -> None:
    """On a consistent model the two routines agree, with nothing blocked."""
    from_me = spectra_me(toy_model, ["R3"], seed=0)
    from_ccme, blocked_core = spectra_ccme(toy_model, ["R3"], seed=0)
    assert blocked_core == []
    assert _ids(from_me) == _ids(from_ccme)


@pytest.mark.parametrize("problem_type", (MIN_NET_MILP, TRADE_OFF))
def test_pathway_exclusion_gives_distinct_models(
    toy_model: Model, problem_type: str
) -> None:
    """Pathway exclusion returns different models, each keeping the core.

    With ``R5`` as the only core reaction the toy model has three ways to
    satisfy it, so three distinct solutions should be reachable.
    """
    models: List[Model] = spectra_me(
        toy_model,
        ["R5"],
        n_solutions=3,
        alt_solution_method=PATHWAY_EXCLUSION,
        problem_type=problem_type,
        seed=0,
    )
    assert len(models) == 3
    reaction_sets = [frozenset(_ids(m)) for m in models]
    assert len(set(reaction_sets)) == 3
    for model in models:
        assert "R5" in _ids(model)


def test_pathway_exclusion_needs_a_milp(toy_model: Model) -> None:
    """Pathway exclusion with an LP formulation is an error."""
    with pytest.raises(SpectraError, match="mixed-integer"):
        spectra_me(
            toy_model,
            ["R5"],
            n_solutions=2,
            alt_solution_method=PATHWAY_EXCLUSION,
            problem_type=MIN_NET_LP,
            seed=0,
        )


def test_core_direction_gives_several_models(toy_model: Model) -> None:
    """The coreDirection method returns as many models as asked for."""
    models = spectra_me(toy_model, ["R3"], n_solutions=3, seed=0)
    assert len(models) == 3
    for model in models:
        assert "R3" in _ids(model)


def test_single_solution_returns_a_bare_model(toy_model: Model) -> None:
    """One solution comes back as a model, not a list of one."""
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    assert isinstance(extracted, Model)


def test_accepts_core_as_indices_and_objects(toy_model: Model) -> None:
    """Core reactions may be given as identifiers, indices or objects."""
    by_id = _ids(spectra_me(toy_model, ["R3"], seed=0))
    by_index = _ids(spectra_me(toy_model, [2], seed=0))
    by_object = _ids(spectra_me(toy_model, [toy_model.reactions.R3], seed=0))
    assert by_id == by_index == by_object


def test_rejects_unknown_core_reaction(toy_model: Model) -> None:
    """A core reaction that is not in the model is an error."""
    with pytest.raises(SpectraError, match="not"):
        spectra_me(toy_model, ["nope"], seed=0)


def test_rejects_unknown_weighted_reaction(toy_model: Model) -> None:
    """A weight for a reaction that is not in the model is an error."""
    with pytest.raises(SpectraError, match="not in the model"):
        spectra_me(toy_model, ["R3"], weights={"nope": 1.0}, seed=0)


def test_topology_mode_extracts_from_an_inconsistent_model(
    blocked_model: Model,
) -> None:
    """Under the accumulation condition the dead-end core works."""
    extracted = spectra_me(blocked_model, ["R3"], consistency_type="topology", seed=0)
    assert {"R1", "R2", "R3"} <= _ids(extracted)


def test_does_not_modify_the_input_model(toy_model: Model) -> None:
    """Extraction leaves the input model's bounds and variables untouched."""
    before_bounds = {r.id: r.bounds for r in toy_model.reactions}
    n_variables = len(toy_model.variables)

    spectra_me(toy_model, ["R3"], problem_type=MIN_NET_MILP, seed=0)

    assert {r.id: r.bounds for r in toy_model.reactions} == before_bounds
    assert len(toy_model.variables) == n_variables


def test_remove_genes_drops_orphans(toy_model: Model) -> None:
    """Orphaned genes are dropped when asked for."""
    toy_model.reactions.R7.gene_reaction_rule = "g_unused"
    toy_model.reactions.R3.gene_reaction_rule = "g_used"

    kept = spectra_me(toy_model, ["R3"], remove_genes=False, seed=0)
    pruned = spectra_me(toy_model, ["R3"], remove_genes=True, seed=0)

    # R7 is not in either extracted model, so its gene is an orphan.
    assert "R7" not in _ids(kept) and "R7" not in _ids(pruned)
    assert "g_used" in {g.id for g in pruned.genes}
    assert len({g.id for g in pruned.genes}) <= len({g.id for g in kept.genes})
