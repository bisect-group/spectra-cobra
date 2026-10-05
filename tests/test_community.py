"""Test joining organism models into a community."""

import pytest
from cobra import Metabolite, Model, Reaction

from spectra_cobra import SpectraError
from spectra_cobra.community import (
    COUPLING_THRESHOLD,
    ORGANISM_SEPARATOR,
    POOLED,
    SHARED,
    build_community_model,
)


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

    assert growth <= COUPLING_THRESHOLD + 1e-9, (
        "with the producer shut down the consumer should be held to the "
        "coupling threshold, which is all a dead unit may carry"
    )


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


def test_shared_mode_needs_no_transports_and_agrees_with_pooled(pair) -> None:
    """Sharing the external compartment is the same model, built smaller."""
    medium = {"glc_e": (-10.0, 1000.0)}
    pooled = build_community_model(list(pair), organisms=["P", "C"], pool_medium=medium)
    shared = build_community_model(
        list(pair), organisms=["P", "C"], mode=SHARED, pool_medium=medium
    )

    assert not [r for r in shared.model.reactions if r.id.startswith("TR_")]
    assert [r for r in pooled.model.reactions if r.id.startswith("TR_")]
    assert len(shared.model.reactions) < len(pooled.model.reactions)

    def consumer_growth(community):
        model = community.model
        model.objective = model.reactions.get_by_id(community.biomass_reactions["C"])
        return model.slim_optimize()

    assert consumer_growth(shared) == pytest.approx(consumer_growth(pooled), rel=1e-6)


def test_shared_mode_allows_only_one_pool(pair) -> None:
    """Its pool is the external compartment, and there is only one of those."""
    with pytest.raises(SpectraError, match="only be one pool"):
        build_community_model(
            list(pair),
            organisms=["P", "C"],
            mode=SHARED,
            pools={"a": ["P"], "b": ["C"]},
            environment=["a"],
        )


def test_shared_mode_rejects_link_bounds(pair) -> None:
    """There is no per-unit reaction left to bound."""
    with pytest.raises(SpectraError, match="no such reactions to bound"):
        build_community_model(
            list(pair),
            organisms=["P", "C"],
            mode=SHARED,
            link_bounds={"glc_e": (0.0, 1.0)},
        )


@pytest.fixture(scope="function")
def three(solver: str):
    """Return three units: two on the blood, one reachable only through B."""
    return (
        _organism("a", solver, eats="glc", makes="ac"),
        _organism("b", solver, eats="glc", makes="ac"),
        _organism("c", solver, eats="ac", makes="co2"),
    )


def _blood_and_interface(three, **kwargs):
    """Join three units with one open pool and one closed interface."""
    return build_community_model(
        list(three),
        organisms=["A", "B", "C"],
        pools={"Bl": ["A", "B"], "iface": ["B", "C"]},
        environment=["Bl"],
        pool_medium={"glc_e": (-10.0, 1000.0)},
        **kwargs,
    )


def test_an_interface_pool_gets_no_exchanges(three) -> None:
    """Only the pools named in environment may reach outside."""
    community = _blood_and_interface(three)

    assert "EX_glc_Bl" in community.model.reactions
    assert not [
        rxn_id for rxn_id in community.community_exchanges if rxn_id.endswith("_iface")
    ]


def test_a_unit_off_the_open_pool_lives_through_its_neighbour(three) -> None:
    """C never touches the blood, so everything it eats comes from B."""
    community = _blood_and_interface(three)
    model = community.model
    model.objective = model.reactions.get_by_id(community.biomass_reactions["C"])

    assert (model.slim_optimize() or 0.0) > 1e-6

    with model:
        model.reactions.get_by_id(community.biomass_reactions["B"]).bounds = (0.0, 0.0)
        starved = model.slim_optimize() or 0.0
    assert starved <= COUPLING_THRESHOLD + 1e-9


def test_several_pools_must_say_which_reach_the_environment(three) -> None:
    """Defaulting to all of them would let an interface bypass the blood."""
    with pytest.raises(SpectraError, match="which of them"):
        build_community_model(
            list(three),
            organisms=["A", "B", "C"],
            pools={"Bl": ["A", "B"], "iface": ["B", "C"]},
        )


@pytest.mark.parametrize(
    "pools, environment, message",
    [
        ({"Bl": ["A", "B", "nope"]}, ["Bl"], "do not exist"),
        ({"Bl": ["A", "B"]}, ["Bl"], "in at least one pool"),
        ({"A": ["A", "B", "C"]}, ["A"], "cannot be named after a unit"),
        ({"Bl": ["A", "B", "C"]}, ["missing"], "unknown pools"),
    ],
)
def test_pools_are_validated(three, pools, environment, message) -> None:
    """A pool graph that does not make sense should say so, not solve."""
    with pytest.raises(SpectraError, match=message):
        build_community_model(
            list(three), organisms=["A", "B", "C"], pools=pools, environment=environment
        )


def test_link_bounds_constrain_what_crosses_into_a_pool(three) -> None:
    """This is where a measured uptake rate goes."""
    community = _blood_and_interface(three, link_bounds={"ac_e": (-2.0, 3.0)})
    model = community.model

    assert model.reactions.get_by_id("TR_ac_Bl__A").bounds == (-2.0, 3.0)
    assert model.reactions.get_by_id("TR_ac_iface__C").bounds == (-2.0, 3.0)
    assert model.reactions.get_by_id("TR_glc_Bl__A").bounds == (-1000.0, 1000.0)


def test_pool_medium_constrains_the_community_exchanges(three) -> None:
    """Anything unnamed stays shut to uptake."""
    community = _blood_and_interface(three)
    model = community.model

    assert model.reactions.EX_glc_Bl.bounds == (-10.0, 1000.0)
    assert model.reactions.EX_ac_Bl.lower_bound == 0.0


def test_a_closed_community_is_refused(solver: str) -> None:
    """With no way in or out nothing can grow, so say so at build time.

    The unit here has no demand or sink of its own, so shutting the pool
    really does close it. A unit that kept one would stay open, which is
    why the check looks for any boundary reaction rather than for an
    exchange on a pool.
    """
    model = Model("sealed")
    outside = Metabolite("glc_e", compartment="e")
    inside = Metabolite("x_c", compartment="c")
    model.add_metabolites([outside, inside])
    uptake = Reaction("UP", lower_bound=0.0, upper_bound=1000.0)
    back = Reaction("BACK", lower_bound=0.0, upper_bound=1000.0)
    model.add_reactions([uptake, back])
    uptake.add_metabolites({outside: -1.0, inside: 1.0})
    back.add_metabolites({inside: -1.0, outside: 1.0})
    model.add_boundary(outside, type="exchange")
    model.objective = uptake
    model.solver = solver

    with pytest.raises(SpectraError, match="no boundary reaction"):
        build_community_model([model], organisms=["X"], environment=[])


def _database(solver: str) -> Model:
    """Return a database holding the producer's reactions and one more."""
    model = _organism("db", solver, eats="glc", makes="ac")
    spare = Reaction("SPARE", lower_bound=0.0, upper_bound=1000.0)
    model.add_reactions([spare])
    spare.add_metabolites({model.metabolites.get_by_id("x_c"): -1.0})
    return model


def test_a_database_is_folded_in_and_its_extra_reactions_recorded(
    pair, solver: str
) -> None:
    """The unit holds the whole database; the draft supplies the bounds."""
    producer, consumer = pair
    producer.reactions.CONV.bounds = (0.0, 7.0)

    community = build_community_model(
        [producer, consumer],
        organisms=["P", "C"],
        databases={"P": _database(solver)},
    )

    assert "SPARE__P" in community.model.reactions
    assert community.database_reactions["P"] == ("SPARE__P",)
    assert "CONV__P" in community.draft_reactions["P"]
    assert community.model.reactions.get_by_id("CONV__P").bounds == (
        0.0,
        7.0,
    ), "where both have the reaction, what is known about this unit wins"
    assert "C" not in community.database_reactions


def test_the_anchor_is_reported_with_how_it_was_decided(pair) -> None:
    """The caller has to be able to see what their model got coupled to."""
    community = build_community_model(list(pair), organisms=["P", "C"])

    assert community.anchor_sources == {"P": "objective", "C": "objective"}
    assert community.anchor_reactions is community.biomass_reactions
    assert "biomass__P" in community.coupling_summary()

    given = build_community_model(
        list(pair),
        organisms=["P", "C"],
        biomass_reactions={"P": "CONV", "C": "biomass"},
    )
    assert given.biomass_reactions["P"] == "CONV__P"
    assert given.anchor_sources["P"] == "given"


def test_an_ambiguous_anchor_is_refused_rather_than_guessed(solver: str) -> None:
    """Coupling a unit to the wrong reaction would not fail loudly."""
    model = Model("ambiguous")
    inside = Metabolite("x_c", compartment="c")
    model.add_metabolites([inside])
    first = Reaction("GROW", lower_bound=0.0, upper_bound=1000.0)
    second = Reaction("ATPM", lower_bound=0.0, upper_bound=1000.0)
    model.add_reactions([first, second])
    first.add_metabolites({inside: -1.0})
    second.add_metabolites({inside: -1.0})
    model.objective = {first: 1.0, second: 1.0}
    model.solver = solver

    with pytest.raises(SpectraError, match="candidates look like"):
        build_community_model([model], organisms=["X"])


def test_decompose_gives_each_unit_back_untagged(pair) -> None:
    """What was a transport into the pool comes back as an exchange."""
    community = build_community_model(
        list(pair), organisms=["P", "C"], pool_medium={"glc_e": (-10.0, 1000.0)}
    )

    parts = community.decompose()

    assert set(parts) == {"P", "C"}
    producer = parts["P"]
    # SEC is written "ac_e ->", which the builder rightly treats as an
    # exchange and replaces with a transport; that transport is what comes
    # back as this unit's own exchange.
    assert {"UP", "CONV", "biomass"} <= {r.id for r in producer.reactions}
    assert not [r.id for r in producer.reactions if ORGANISM_SEPARATOR in r.id]
    assert not [m.id for m in producer.metabolites if ORGANISM_SEPARATOR in m.id]
    assert producer.boundary, "the unit must still be able to feed itself"


def test_multi_tissue_replicates_one_model_across_tissues(solver: str) -> None:
    """Tissues differ by their bounds and their pool, not by their network."""
    from spectra_cobra import build_multi_tissue_model

    gem = _organism("gem", solver, eats="glc", makes="ac")

    body = build_multi_tissue_model(
        gem,
        tissues=["tis1", "tis2", "tis3"],
        pools={"Bl": ["tis1", "tis2"], "tis2_tis3": ["tis2", "tis3"]},
        environment=["Bl"],
        anchor_reactions={t: "biomass" for t in ("tis1", "tis2", "tis3")},
        pool_medium={"glc_e": (-10.0, 1000.0)},
    )

    assert body.organisms == ("tis1", "tis2", "tis3")
    assert all(source == "given" for source in body.anchor_sources.values())
    assert "EX_glc_Bl" in body.model.reactions
    assert "EX_glc_tis2_tis3" not in body.model.reactions
    # tis3 never touches the blood, so it eats only what tis2 passes on.
    assert not [
        rxn_id for rxn_id in body.reactions_of["tis3"] if rxn_id.endswith("_Bl__tis3")
    ]


def _eater(name: str, met_id: str, solver: str, compartment: str = "e") -> Model:
    """Return a one-reaction unit that eats one external metabolite."""
    model = Model(name)
    outside = Metabolite(met_id, compartment=compartment)
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


@pytest.mark.parametrize("mode", [POOLED, SHARED])
def test_both_spellings_meet_in_the_pool_in_either_mode(solver: str, mode) -> None:
    """ "glc_e" and "glc[e]" are the same compound and must share a row.

    Getting this wrong is silent and total: each unit would trade through
    a pool of its own, named after its own spelling, and the community
    would never cross-feed.
    """
    community = build_community_model(
        [_eater("a", "glc_e", solver), _eater("b", "glc[e]", solver)],
        organisms=["A", "B"],
        mode=mode,
        check=False,
    )

    pooled = [m.id for m in community.model.metabolites if m.compartment in {"u", "e"}]
    assert len(pooled) == 1, f"got {pooled}"
    assert len(community.community_exchanges) == 1


def test_the_external_compartment_is_found_per_unit(solver: str) -> None:
    """The units need not agree on what it is called."""
    community = build_community_model(
        [
            _eater("a", "glc_e", solver),
            _eater("b", "glc_extracellular", solver, compartment="extracellular"),
        ],
        organisms=["A", "B"],
        check=False,
    )

    assert [m.id for m in community.model.metabolites if m.compartment == "u"] == [
        "glc_u"
    ], "both units' glucose should reach the same pool metabolite"


def test_units_sharing_nothing_are_reported(solver: str, caplog) -> None:
    """Identifiers from different namespaces cannot be guessed apart."""
    with caplog.at_level("WARNING"):
        build_community_model(
            [_eater("a", "glc__D_e", solver), _eater("b", "glc_D[e]", solver)],
            organisms=["A", "B"],
            check=False,
        )

    assert "nothing can cross-feed" in caplog.text
    assert "metabolite_key" in caplog.text


def test_metabolite_key_reconciles_two_namespaces(solver: str) -> None:
    """The escape hatch for models that do not share a naming scheme."""
    import re

    community = build_community_model(
        [_eater("a", "glc__D_e", solver), _eater("b", "glc_D[e]", solver)],
        organisms=["A", "B"],
        metabolite_key=lambda met: re.sub(r"_+D(_e|\[e\])$", "_D", met.id),
        check=False,
    )

    assert [m.id for m in community.model.metabolites if m.compartment == "u"] == [
        "glc_D_u"
    ]


def test_a_medium_naming_an_absent_metabolite_says_so(solver: str, caplog) -> None:
    """A typo in a medium should not quietly do nothing."""
    with caplog.at_level("WARNING"):
        build_community_model(
            [_eater("a", "glc_e", solver)],
            organisms=["A"],
            pool_medium={"glc_e": (-10.0, 1000.0), "nosuch_e": (-5.0, 0.0)},
            check=False,
        )

    assert "are in no unit and were ignored" in caplog.text
    assert "nosuch_e" in caplog.text


def test_a_unit_that_merely_lost_a_tie_is_not_called_inert(solver: str, caplog) -> None:
    """Maximising every anchor at once lands on one point of a flat face.

    Two units living on the same nutrient can share it in any proportion,
    so the solver is free to give one of them all of it and the other
    none. Reporting the loser as dead would be a false alarm, and the
    reading is only a lower bound, so a zero is always checked again on
    its own.
    """
    with caplog.at_level("WARNING"):
        community = build_community_model(
            [_eater("a", "glc_e", solver), _eater("b", "glc_e", solver)],
            organisms=["A", "B"],
            pool_medium={"glc_e": (-10.0, 1000.0)},
            couple=False,
        )

    assert "inert" not in caplog.text, caplog.text
    assert min(community.anchor_capacity.values()) > 0.0


def test_a_unit_that_really_is_inert_is_reported(solver: str, caplog) -> None:
    """The contrast: nothing on offer feeds this one at all."""
    with caplog.at_level("WARNING"):
        community = build_community_model(
            [_eater("a", "glc_e", solver), _eater("b", "xyz_e", solver)],
            organisms=["A", "B"],
            pool_medium={"glc_e": (-10.0, 1000.0)},
            couple=False,
        )

    assert "is inert" in caplog.text
    assert community.anchor_capacity["B"] == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("problem_type", ["minNetLP", "minNetMILP"])
def test_a_multi_tissue_model_can_be_extracted_and_split_again(
    solver: str, problem_type
) -> None:
    """The body is a universal model like any other, in either formulation."""
    from spectra_cobra import build_multi_tissue_model, spectra_me

    gem = _organism("gem", solver, eats="glc", makes="ac")
    tissues = ["tis1", "tis2", "tis3"]
    body = build_multi_tissue_model(
        gem,
        tissues=tissues,
        pools={"Bl": ["tis1", "tis2"], "tis2_tis3": ["tis2", "tis3"]},
        environment=["Bl"],
        anchor_reactions={t: "biomass" for t in tissues},
        pool_medium={"glc_e": (-10.0, 1000.0)},
        couple=False,
    )

    core = [f"biomass{ORGANISM_SEPARATOR}{t}" for t in tissues]
    extracted = spectra_me(
        body.model, core, tol=1e-4, problem_type=problem_type, seed=0
    )

    assert all(rxn_id in extracted.reactions for rxn_id in core)
    assert len(extracted.reactions) < len(body.model.reactions)

    parts = body.with_model(extracted).decompose()
    assert set(parts) == set(tissues)
    for tissue, part in parts.items():
        assert "biomass" in part.reactions, tissue
        assert not [r.id for r in part.reactions if ORGANISM_SEPARATOR in r.id]
