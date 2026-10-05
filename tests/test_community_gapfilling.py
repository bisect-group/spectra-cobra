"""Test gap-filling every unit of a community at once.

The toy is a four-step chain, ``glc -> a -> b -> c``, ending in an anchor
that consumes ``c``. Two drafts each hold part of the chain:

* ``A`` can take up glucose and make ``b``, and can export ``b``, but
  cannot turn ``b`` into ``c``.
* ``B`` can turn ``b`` into ``c``, but cannot take up glucose or make
  ``b`` at all.

Separately, ``A`` must borrow one reaction and ``B`` must borrow two.
Together, ``A`` still needs its one -- it has to feed itself -- but ``B``
needs none, because what it is missing arrives from ``A`` through the
pool. That difference is the whole argument for gap-filling a community
rather than its members one at a time.
"""

import pytest
from cobra import Metabolite, Model, Reaction

from spectra_cobra import (
    MIN_NET_LP,
    MIN_NET_MILP,
    SpectraError,
    build_community_model,
    spectra_me,
)
from spectra_cobra.gapfilling import (
    LP_EXCHANGE_WEIGHT,
    MILP_EXCHANGE_WEIGHT,
    _exchange_weight,
    gapfill_community,
)

#: Every reaction of the chain, as (id, {metabolite: coefficient}).
CHAIN = {
    "UP": {"glc_e": -1.0, "a_c": 1.0},
    "R1": {"a_c": -1.0, "b_c": 1.0},
    "EXP": {"b_c": -1.0, "b_e": 1.0},
    "R2": {"b_c": -1.0, "c_c": 1.0},
    "BIO": {"c_c": -1.0},
}

#: What each draft already has. The rest it must borrow or be given.
DRAFTS = {
    "A": ("UP", "R1", "EXP", "BIO"),
    "B": ("EXP", "R2", "BIO"),
}


def _build(name: str, reactions, solver: str) -> Model:
    """Return a model holding the named reactions of the chain.

    Parameters
    ----------
    name : str
        The model's identifier.
    reactions : iterable of str
        Which reactions of `CHAIN` to include.
    solver : str
        The solver to attach.

    Returns
    -------
    cobra.Model
        The model, with exchanges on its external metabolites and ``BIO``
        as its objective.

    """
    model = Model(name)
    wanted = list(reactions)
    met_ids = sorted({m for r in wanted for m in CHAIN[r]})
    model.add_metabolites(
        [Metabolite(m, compartment=m.rsplit("_", 1)[-1]) for m in met_ids]
    )
    for rxn_id in wanted:
        reaction = Reaction(
            rxn_id,
            lower_bound=-1000.0 if rxn_id == "EXP" else 0.0,
            upper_bound=1000.0,
        )
        model.add_reactions([reaction])
        reaction.add_metabolites(
            {model.metabolites.get_by_id(m): c for m, c in CHAIN[rxn_id].items()}
        )
    for met in model.metabolites:
        if met.compartment == "e":
            model.add_boundary(met, type="exchange")
    model.objective = model.reactions.BIO
    model.solver = solver
    return model


@pytest.fixture(scope="function")
def database(solver: str) -> Model:
    """Return the whole chain, which either draft may borrow from."""
    return _build("database", CHAIN, solver)


@pytest.fixture(scope="function")
def community(solver: str, database: Model):
    """Return the two drafts joined, each backed by the whole chain."""
    drafts = [_build(name, ids, solver) for name, ids in DRAFTS.items()]
    return build_community_model(
        drafts,
        organisms=list(DRAFTS),
        databases={name: database.copy() for name in DRAFTS},
        pool_medium={"glc_e": (-10.0, 1000.0)},
    )


def test_the_database_shows_up_as_candidates_not_as_the_draft(community) -> None:
    """Only what a unit lacks is a candidate to be added."""
    assert {r.rsplit("__", 1)[0] for r in community.database_reactions["A"]} == {"R2"}
    assert {r.rsplit("__", 1)[0] for r in community.database_reactions["B"]} == {
        "UP",
        "R1",
    }


@pytest.mark.parametrize("problem_type", [MIN_NET_MILP, MIN_NET_LP])
def test_a_unit_is_fed_rather_than_gap_filled(community, problem_type) -> None:
    """B needs nothing, because A's secretion is what it was missing."""
    result = gapfill_community(community, problem_type=problem_type, seed=0)

    assert result.added["A"] == ("R2",), "A has to feed itself whatever B does"
    assert result.added["B"] == (), "B lives off the pool instead of borrowing"


def test_every_unit_grows_in_the_gap_filled_community(community) -> None:
    """Making each anchor core is what stops a unit being left for dead."""
    result = gapfill_community(community, seed=0)
    model = result.community.model

    for unit in result.community.organisms:
        anchor = model.reactions.get_by_id(result.community.biomass_reactions[unit])
        with model:
            model.objective = anchor
            assert (model.slim_optimize() or 0.0) > 1e-6, f"{unit} cannot grow"


def test_filling_the_units_one_at_a_time_costs_more(community, database) -> None:
    """The comparison the community formulation exists to win."""
    together = gapfill_community(community, seed=0)
    jointly = sum(len(ids) for ids in together.added.values())

    alone = 0
    for name, draft_ids in DRAFTS.items():
        universal = database.copy()
        for exchange in universal.boundary:
            exchange.lower_bound = -10.0 if "glc" in exchange.id else 0.0
        weights = {
            rxn.id: (0.0 if rxn.id in set(draft_ids) else 1.0)
            for rxn in universal.reactions
        }
        filled = spectra_me(
            universal, ["BIO"], tol=1e-4, weights=weights, problem_type=MIN_NET_MILP
        )
        alone += len({r.id for r in filled.reactions} - set(draft_ids))

    assert jointly < alone, f"community added {jointly}, one at a time added {alone}"


def test_decomposing_gives_back_models_that_grow(community) -> None:
    """A gap-filled unit has to stand on its own once it is handed back."""
    result = gapfill_community(community, seed=0)

    models = result.models()

    assert set(models) == {"A", "B"}
    for name, model in models.items():
        assert not [r.id for r in model.reactions if "__" in r.id]
        assert "BIO" in model.reactions
        # Standing alone it may need the pool opened, which decompose
        # turns back into this unit's own exchange reactions.
        for exchange in model.boundary:
            exchange.lower_bound = min(exchange.lower_bound, -10.0)
        assert (model.slim_optimize() or 0.0) > 1e-6, f"{name} cannot grow alone"


def test_an_unknown_core_reaction_is_refused(community) -> None:
    """A typo should not quietly drop the requirement."""
    with pytest.raises(SpectraError, match="not in the community model"):
        gapfill_community(community, core_reactions=["NOPE"])


def test_without_core_anchors_a_unit_can_be_left_out(community) -> None:
    """The contrast that shows what anchors_are_core is buying."""
    result = gapfill_community(community, anchors_are_core=False, seed=0)

    assert not result.core
    assert (
        sum(len(ids) for ids in result.added.values()) == 0
    ), "with nothing required, the cheapest community adds nothing at all"


def test_keep_draft_retains_what_the_solution_does_not_use(community) -> None:
    """Gap-filling adds; it must not quietly subtract."""
    kept = gapfill_community(community, seed=0)
    lean = gapfill_community(community, keep_draft=False, seed=0)

    drafted = {r for ids in community.draft_reactions.values() for r in ids}
    present = {r.id for r in kept.community.model.reactions}
    assert drafted <= present, "every draft reaction should survive"

    trimmed = {r.id for r in lean.community.model.reactions}
    assert len(trimmed) <= len(present)
    assert lean.added == kept.added, "what was borrowed does not change"


@pytest.mark.parametrize(
    "problem_type, expected",
    [(MIN_NET_MILP, MILP_EXCHANGE_WEIGHT), (MIN_NET_LP, LP_EXCHANGE_WEIGHT)],
)
def test_the_exchange_weight_follows_the_formulation(problem_type, expected) -> None:
    """An LP sums flux, so it must not be charged per unit sharing a pool.

    Measured on two hCom organisms: weighting community exchanges 1 under
    ``minNetLP`` adds 61 reactions where weighting them 0 adds 8, which is
    also what the mixed-integer answer comes to.
    """
    assert _exchange_weight(None, problem_type) == expected
    assert _exchange_weight(0.5, problem_type) == 0.5
