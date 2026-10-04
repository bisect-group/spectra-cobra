"""Test gap-filling towards growth and towards metabolic tasks."""

import pytest
from cobra import Model

from spectra_cobra import (
    DEFAULT_MIN_GROWTH,
    MetabolicTask,
    SpectraError,
    gapfill_for_growth,
    gapfill_for_tasks,
)

#: Everything the textbook model needs besides a carbon source.
SALTS = {
    "EX_pi_e": -1000.0,
    "EX_h2o_e": -1000.0,
    "EX_h_e": -1000.0,
    "EX_nh4_e": -1000.0,
    "EX_o2_e": -1000.0,
    "EX_co2_e": -1000.0,
}
#: The reactions knocked out to make a draft with gaps.
GAPS = ["ENO", "ACONTb", "SUCDi", "NH4t"]


@pytest.fixture(scope="function")
def universal(solver: str) -> Model:
    """Return the textbook model, which stands in for a universal model."""
    from cobra.io import load_model

    model = load_model("textbook")
    model.solver = solver
    return model


@pytest.fixture(scope="function")
def draft(universal: Model) -> Model:
    """Return the textbook model with four reactions removed."""
    model = universal.copy()
    model.remove_reactions(GAPS, remove_orphans=False)
    return model


@pytest.fixture(scope="function")
def media() -> dict:
    """Return six single-carbon-source media."""
    return {
        carbon: {**SALTS, f"EX_{carbon}_e": -10.0}
        for carbon in ("glc__D", "ac", "succ", "mal__L", "pyr", "akg")
    }


def test_the_default_growth_demand_is_a_real_rate() -> None:
    """A numerical floor would be satisfied by a model that barely lives."""
    assert DEFAULT_MIN_GROWTH == 0.1


def test_gapfilling_restores_growth_in_every_medium(draft, universal, media) -> None:
    """The draft grows nowhere; the gap-filled model grows everywhere."""
    result = gapfill_for_growth(draft, universal, media, seed=0)

    assert set(result.satisfied) == set(media)
    assert result.unsatisfied == ()
    assert not result.failed
    assert set(result.added) <= set(GAPS)


def test_the_draft_is_not_modified(draft, universal, media) -> None:
    """Gap-filling returns a new model rather than editing the input."""
    before = {rxn.id for rxn in draft.reactions}

    gapfill_for_growth(draft, universal, media, seed=0)

    assert {rxn.id for rxn in draft.reactions} == before


def test_growth_is_demanded_through_a_bound_not_through_tol(
    draft, universal, media
) -> None:
    """The result must not depend on tol, which is what the bound buys.

    Setting tol to the growth rate instead would also demand that rate of
    every other core reaction, and would make the answer move with tol.
    """
    one = gapfill_for_growth(draft, universal, media, tol=1e-4, seed=0)
    two = gapfill_for_growth(draft, universal, media, tol=1e-6, seed=0)

    assert set(one.added) == set(two.added)


def test_a_higher_growth_demand_is_honoured(draft, universal, media) -> None:
    """Asking for more growth has to produce a model that delivers it."""
    glucose = {"glc__D": media["glc__D"]}

    result = gapfill_for_growth(draft, universal, glucose, min_growth=0.5, seed=0)

    assert result.satisfied == ("glc__D",)
    model = result.model
    with model:
        for reaction in model.boundary:
            reaction.lower_bound = 0.0
        for rxn_id, lower in glucose["glc__D"].items():
            model.reactions.get_by_id(rxn_id).lower_bound = lower
        assert model.slim_optimize() >= 0.5


def test_min_growth_none_only_asks_for_some_growth(draft, universal, media) -> None:
    """Opting out of the bound still requires the model to be viable."""
    result = gapfill_for_growth(
        draft, universal, {"glc__D": media["glc__D"]}, min_growth=None, seed=0
    )

    assert result.satisfied == ("glc__D",)


def test_trade_off_is_refused_for_growth(draft, universal, media) -> None:
    """It requires every included reaction to carry tol, which growth cannot."""
    with pytest.raises(SpectraError, match="tradeOff requires"):
        gapfill_for_growth(draft, universal, media, problem_type="tradeOff")


def test_an_impossible_medium_is_reported_not_raised(draft, universal, media) -> None:
    """One unusable condition must not lose the others."""
    impossible = dict(media)
    # No carbon source at all, so nothing can grow on it.
    impossible["empty"] = dict(SALTS)

    result = gapfill_for_growth(draft, universal, impossible, seed=0)

    assert "empty" in result.unsatisfied or "empty" in result.failed
    assert "glc__D" in result.satisfied


def test_a_medium_naming_an_unknown_reaction_raises(draft, universal) -> None:
    """A typo in a medium should be reported immediately."""
    with pytest.raises(SpectraError, match="does not have"):
        gapfill_for_growth(draft, universal, [{"EX_nope_e": -10.0}])


def test_free_reactions_cost_nothing_to_include(draft, universal, media) -> None:
    """Naming a reaction free should never make the answer worse."""
    result = gapfill_for_growth(draft, universal, media, free_reactions=GAPS, seed=0)

    assert set(result.satisfied) == set(media)


def test_gapfilling_towards_tasks(universal) -> None:
    """A draft that fails its tasks is repaired until it does not."""
    from spectra_cobra import TaskEquation

    atp = TaskEquation(
        {"atp_c": -1, "h2o_c": -1, "adp_c": 1, "pi_c": 1, "h_c": 1}, (1.0, 1000.0)
    )
    medium = {
        "glc__D_e": (0, 1000),
        "o2_e": (0, 1000),
        "nh4_e": (0, 1000),
        "pi_e": (0, 1000),
        "h2o_e": (0, 1000),
        "h_e": (0, 1000),
    }
    waste = {"co2_e": (0, 1000), "h2o_e": (0, 1000), "h_e": (0, 1000)}
    tasks = [
        MetabolicTask("atp", inputs=medium, outputs=waste, equations=(atp,)),
        MetabolicTask(
            "glutamate", inputs=medium, outputs={"glu__L_c": (0.1, 1000), **waste}
        ),
    ]
    draft = universal.copy()
    draft.remove_reactions(["CYTBD", "GLUSy", "AKGDH"], remove_orphans=False)

    result = gapfill_for_tasks(draft, universal, tasks, seed=0)

    assert result.unsatisfied == (), result.summary()
    assert set(result.satisfied) == {"atp", "glutamate"}
