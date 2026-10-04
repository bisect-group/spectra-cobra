"""Test reducing a community to a minimal microbiome."""

import pytest
from cobra import Metabolite, Model, Reaction

from spectra_cobra import SpectraError, build_community_model, minimal_microbiome


def _organism(name: str, eats: str, makes: str) -> Model:
    """Return a one-step organism: eats one metabolite, secretes another.

    Parameters
    ----------
    name : str
        The model's identifier.
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
    biomass.add_metabolites({inside: -1.0})
    secrete.add_metabolites({waste: -1.0})
    model.add_boundary(food, type="exchange")
    model.add_boundary(waste, type="exchange")
    model.objective = biomass
    return model


@pytest.fixture(scope="function")
def gut(solver: str):
    """Return a three-member community on glucose.

    ``P`` turns glucose into acetate, ``C`` turns acetate into CO2, and
    ``F`` turns glucose into lactate without helping anyone. So CO2
    production needs both P and C, lactate needs only F, and growth alone
    needs just one of them.
    """
    community = build_community_model(
        [
            _organism("producer", "glc", "ac"),
            _organism("consumer", "ac", "co2"),
            _organism("freeloader", "glc", "lac"),
        ],
        organisms=["P", "C", "F"],
    )
    model = community.model
    model.solver = solver
    for exchange in community.community_exchanges:
        model.reactions.get_by_id(exchange).lower_bound = 0.0
    model.reactions.EX_glc_u.lower_bound = -10.0
    return community


def test_the_product_decides_who_is_needed(gut) -> None:
    """Requiring CO2 keeps the chain that makes it, and nothing else."""
    result = minimal_microbiome(gut, products=["EX_co2_u"], seed=0)

    assert set(result.present) == {"P", "C"}
    assert result.absent == ("F",)


def test_a_different_product_keeps_a_different_member(gut) -> None:
    """Lactate comes from the freeloader, so now it is the one that matters."""
    result = minimal_microbiome(gut, products=["EX_lac_u"], seed=0)

    assert set(result.present) == {"F"}


def test_growth_alone_keeps_the_fewest(gut) -> None:
    """With no product to preserve, one organism is enough."""
    result = minimal_microbiome(gut, seed=0)

    assert len(result.present) == 1


def test_the_requirements_are_actually_met(gut) -> None:
    """A smaller community is only correct if it still does the job."""
    result = minimal_microbiome(gut, products=["EX_co2_u"], seed=0)

    assert result.growth >= result.required_growth - 1e-6
    assert result.production >= result.required_production - 1e-6


def test_required_organisms_are_kept(gut) -> None:
    """An organism the caller insists on stays, however useless."""
    result = minimal_microbiome(gut, products=["EX_co2_u"], required=["F"], seed=0)

    assert "F" in result.present
    assert result.membership["F"]


def test_membership_covers_every_organism(gut) -> None:
    """The vector is over the original community, not the surviving one."""
    result = minimal_microbiome(gut, products=["EX_co2_u"], seed=0)

    assert set(result.membership) == {"P", "C", "F"}
    assert [o for o, kept in result.membership.items() if kept] == list(result.present)


def test_the_returned_community_drops_the_absent(gut) -> None:
    """The minimal model must not still contain the organisms it excluded."""
    result = minimal_microbiome(gut, products=["EX_co2_u"], seed=0)

    assert set(result.community.organisms) == set(result.present)
    assert not [
        rxn.id for rxn in result.community.model.reactions if rxn.id.endswith("__F")
    ]


def test_a_stricter_product_demand_keeps_more(gut) -> None:
    """Asking for all of the production cannot keep fewer organisms."""
    lenient = minimal_microbiome(
        gut, products=["EX_co2_u"], product_fraction=0.1, seed=0
    )
    strict = minimal_microbiome(
        gut, products=["EX_co2_u"], product_fraction=1.0, seed=0
    )

    assert len(strict.present) >= len(lenient.present)


def test_an_unknown_product_is_rejected(gut) -> None:
    """A typo should not silently become "no product required"."""
    with pytest.raises(SpectraError, match="No reaction"):
        minimal_microbiome(gut, products=["EX_nope_u"])


def test_an_unknown_required_organism_is_rejected(gut) -> None:
    """Likewise for an organism the community does not have."""
    with pytest.raises(SpectraError, match="Not organisms"):
        minimal_microbiome(gut, required=["nobody"])


def test_a_community_that_cannot_grow_is_reported(gut) -> None:
    """There is no minimal subset of a community that does nothing."""
    gut.model.reactions.EX_glc_u.lower_bound = 0.0

    with pytest.raises(SpectraError, match="cannot grow"):
        minimal_microbiome(gut, products=["EX_co2_u"])


def test_the_input_community_is_not_modified(gut) -> None:
    """Reducing a community must leave the original intact."""
    before = {rxn.id for rxn in gut.model.reactions}
    bounds = {rxn.id: rxn.bounds for rxn in gut.model.reactions}

    minimal_microbiome(gut, products=["EX_co2_u"], seed=0)

    assert {rxn.id for rxn in gut.model.reactions} == before
    assert {rxn.id: rxn.bounds for rxn in gut.model.reactions} == bounds
