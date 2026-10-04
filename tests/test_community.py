"""Test joining organism models into a community."""

import pytest
from cobra import Metabolite, Model, Reaction

from spectra_cobra import SpectraError
from spectra_cobra.community import ORGANISM_SEPARATOR, build_community_model


def _organism(name: str, solver: str, eats: str, makes: str) -> Model:
    """Return a one-step organism: eats one metabolite, secretes another.

    Parameters
    ----------
    name : str
        The model's identifier.
    solver : str
        The solver to attach.
    eats : str
        The extracellular metabolite taken up.
    makes : str
        The extracellular metabolite secreted.

    Returns
    -------
    cobra.Model
        The toy organism.

    """
    model = Model(name)
    food = Metabolite(f"{eats}_e", compartment="e")
    waste = Metabolite(f"{makes}_e", compartment="e")
    inside = Metabolite("x_c", compartment="c")
    model.add_metabolites([food, waste, inside])

    uptake = Reaction("UP", lower_bound=0.0, upper_bound=1000.0)
    convert = Reaction("CONV", lower_bound=0.0, upper_bound=1000.0)
    secrete = Reaction("SEC", lower_bound=0.0, upper_bound=1000.0)
    biomass = Reaction("biomass", lower_bound=0.0, upper_bound=1000.0)
    model.add_reactions([uptake, convert, secrete, biomass])
    uptake.add_metabolites({food: -1.0, inside: 1.0})
    convert.add_metabolites({inside: -1.0, waste: 1.0})
    # Growth costs a little of what the organism took in. This consumes
    # without producing, so cobrapy counts it as a boundary reaction, which
    # is exactly the case the builder must not mistake for an exchange.
    biomass.add_metabolites({inside: -1.0})
    secrete.add_metabolites({waste: -1.0})

    # Exchanges, which the community builder replaces with transports.
    model.add_boundary(food, type="exchange")
    model.add_boundary(waste, type="exchange")
    model.objective = biomass
    model.solver = solver
    return model


@pytest.fixture(scope="function")
def pair(solver: str):
    """Return two organisms where the second lives off the first's waste."""
    producer = _organism("producer", solver, eats="glc", makes="ac")
    consumer = _organism("consumer", solver, eats="ac", makes="co2")
    return producer, consumer


def test_joining_keeps_every_organism_separate(pair) -> None:
    """Each organism's reactions are tagged and its own."""
    community = build_community_model(list(pair), organisms=["P", "C"])

    assert community.organisms == ("P", "C")
    assert all(
        rxn_id.endswith(f"{ORGANISM_SEPARATOR}P")
        for rxn_id in community.reactions_of["P"]
    )
    assert not set(community.reactions_of["P"]) & set(community.reactions_of["C"])
    assert community.organism_of(community.biomass_reactions["C"]) == "C"


def test_the_organisms_own_exchanges_are_replaced(pair) -> None:
    """Keeping them would let every organism draw on the environment."""
    community = build_community_model(list(pair), organisms=["P", "C"])

    assert not [
        rxn.id
        for rxn in community.model.reactions
        if rxn.id.startswith("EX_") and ORGANISM_SEPARATOR in rxn.id
    ]
    assert community.community_exchanges


def test_a_closed_community_cannot_grow(pair) -> None:
    """With the environment shut there is nothing to eat."""
    community = build_community_model(list(pair), organisms=["P", "C"])
    for exchange in community.community_exchanges:
        community.model.reactions.get_by_id(exchange).lower_bound = 0.0

    assert (community.model.slim_optimize() or 0.0) == pytest.approx(0.0, abs=1e-6)


def test_cross_feeding_through_the_shared_pool(pair) -> None:
    """The consumer grows on glucose it cannot eat, via the producer."""
    community = build_community_model(list(pair), organisms=["P", "C"])
    model = community.model
    for exchange in community.community_exchanges:
        model.reactions.get_by_id(exchange).lower_bound = 0.0
    model.reactions.EX_glc_u.lower_bound = -10.0

    model.objective = model.reactions.get_by_id(community.biomass_reactions["C"])
    growth = model.slim_optimize()

    assert growth > 1e-6, "the consumer should live off the producer's acetate"


def test_coupling_shuts_an_organism_down_with_its_biomass(pair) -> None:
    """No growth means no flux anywhere in that organism."""
    community = build_community_model(list(pair), organisms=["P", "C"])
    model = community.model
    for exchange in community.community_exchanges:
        model.reactions.get_by_id(exchange).lower_bound = 0.0
    model.reactions.EX_glc_u.lower_bound = -10.0

    producer_biomass = model.reactions.get_by_id(community.biomass_reactions["P"])
    with model:
        producer_biomass.bounds = (0.0, 0.0)
        model.objective = model.reactions.get_by_id(community.biomass_reactions["C"])
        growth = model.slim_optimize() or 0.0

    assert growth == pytest.approx(
        0.0, abs=1e-2
    ), "with the producer shut down nothing should reach the consumer"


def test_without_coupling_a_dead_organism_still_feeds_the_others(pair) -> None:
    """The contrast that makes the coupling worth its constraints."""
    community = build_community_model(list(pair), organisms=["P", "C"], couple=False)
    model = community.model
    for exchange in community.community_exchanges:
        model.reactions.get_by_id(exchange).lower_bound = 0.0
    model.reactions.EX_glc_u.lower_bound = -10.0

    with model:
        model.reactions.get_by_id(community.biomass_reactions["P"]).bounds = (0.0, 0.0)
        model.objective = model.reactions.get_by_id(community.biomass_reactions["C"])
        growth = model.slim_optimize() or 0.0

    assert growth > 1e-6, "uncoupled, a non-growing producer still works for free"


def test_subset_drops_an_organism_and_keeps_the_pool(pair) -> None:
    """The survivors still trade and still face the same environment."""
    community = build_community_model(list(pair), organisms=["P", "C"])

    only_p = community.subset(["P"])

    assert only_p.organisms == ("P",)
    assert "C" not in only_p.biomass_reactions
    assert not [
        r.id for r in only_p.model.reactions if r.id.endswith(f"{ORGANISM_SEPARATOR}C")
    ]
    assert "EX_glc_u" in only_p.model.reactions


def test_subset_rejects_an_unknown_organism(pair) -> None:
    """A typo should not silently return the whole community."""
    community = build_community_model(list(pair), organisms=["P", "C"])

    with pytest.raises(SpectraError, match="Not organisms"):
        community.subset(["P", "nope"])


@pytest.mark.parametrize(
    "names, message",
    [
        (["P", "P"], "unique"),
        (["P", "has__sep"], "must not contain"),
        (["P"], "organism names"),
    ],
)
def test_organism_names_are_validated(pair, names, message) -> None:
    """The names become identifier suffixes, so they have to be usable."""
    with pytest.raises(SpectraError, match=message):
        build_community_model(list(pair), organisms=names)


@pytest.mark.parametrize("spelling", ["{}_e", "{}[e]"])
def test_both_compartment_spellings_share_a_pool(solver: str, spelling) -> None:
    """ "glc_e" and "glc[e]" must both reach the same shared metabolite.

    Model collections differ on this, and getting it wrong does not fail
    loudly: each organism would simply trade through its own private pool
    named after its own spelling, and the community would never cross-feed.
    """
    from cobra import Metabolite, Model, Reaction

    def organism(name: str) -> Model:
        model = Model(name)
        outside = Metabolite(spelling.format("glc"), compartment="e")
        inside = Metabolite("x_c", compartment="c")
        model.add_metabolites([outside, inside])
        uptake = Reaction("UP", lower_bound=0.0, upper_bound=1000.0)
        biomass = Reaction("biomass", lower_bound=0.0, upper_bound=1000.0)
        model.add_reactions([uptake, biomass])
        uptake.add_metabolites({outside: -1.0, inside: 1.0})
        biomass.add_metabolites({inside: -1.0})
        model.add_boundary(outside, type="exchange")
        model.objective = biomass
        model.solver = solver
        return model

    community = build_community_model(
        [organism("a"), organism("b")], organisms=["A", "B"]
    )

    pooled = [m.id for m in community.model.metabolites if m.compartment == "u"]
    assert pooled == ["glc_u"], f"got {pooled}"
    assert "EX_glc_u" in community.model.reactions
