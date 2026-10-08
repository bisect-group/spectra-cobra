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
    """Orphaned genes are dropped when asked for, and kept when not.

    ``remove_orphans`` in cobrapy covers genes and metabolites together, so
    passing it would drop the genes whatever ``remove_genes`` said. The two
    are handled separately for that reason, and this pins the distinction:
    an earlier version of this test compared the counts with ``<=`` and so
    passed while they were always equal.
    """
    toy_model.reactions.R7.gene_reaction_rule = "g_unused"
    toy_model.reactions.R3.gene_reaction_rule = "g_used"

    kept = spectra_me(toy_model, ["R3"], remove_genes=False, seed=0)
    pruned = spectra_me(toy_model, ["R3"], remove_genes=True, seed=0)

    # R7 is in neither extracted model, so its gene is an orphan in both.
    assert "R7" not in _ids(kept) and "R7" not in _ids(pruned)
    assert {g.id for g in pruned.genes} == {"g_used"}
    assert {g.id for g in kept.genes} == {"g_used", "g_unused"}


def test_genes_are_removed_by_default(toy_model: Model) -> None:
    """The default is to prune, so an extraction carries no orphaned genes."""
    toy_model.reactions.R7.gene_reaction_rule = "g_unused"
    toy_model.reactions.R3.gene_reaction_rule = "g_used"

    extracted = spectra_me(toy_model, ["R3"], seed=0)

    assert [g.id for g in extracted.genes if not g.reactions] == []
    assert {g.id for g in extracted.genes} == {"g_used"}


def test_orphaned_metabolites_go_whatever_remove_genes_says(
    toy_model: Model,
) -> None:
    """Keeping the genes must not also keep metabolites nothing touches."""
    for flag in (False, True):
        extracted = spectra_me(toy_model, ["R3"], remove_genes=flag, seed=0)
        orphans = [met.id for met in extracted.metabolites if not met.reactions]
        assert orphans == [], f"remove_genes={flag} left {orphans}"


def test_gene_reaction_rules_survive_pruning(toy_model: Model) -> None:
    """A kept reaction keeps the genes its rule names."""
    toy_model.reactions.R3.gene_reaction_rule = "g_a or g_b"

    extracted = spectra_me(toy_model, ["R3"], remove_genes=True, seed=0)

    assert {g.id for g in extracted.reactions.R3.genes} == {"g_a", "g_b"}
    assert extracted.reactions.R3.gene_reaction_rule == "g_a or g_b"


def test_warns_when_the_solver_tolerance_is_loose(toy_model: Model, caplog) -> None:
    """A solver tolerance close to tol is called out.

    This is the setting that left core reactions present in a Recon3D
    extraction but unable to carry flux, and it is not one a caller would
    think to check, so it warns rather than failing quietly.
    """
    import logging

    toy_model.tolerance = 1e-7
    with caplog.at_level(logging.WARNING, logger="spectra_cobra.extraction"):
        spectra_me(toy_model, ["R3"], tol=1e-4, seed=0)

    assert any("model.tolerance" in r.message for r in caplog.records)


def test_no_warning_when_the_tolerance_is_tight(toy_model: Model, caplog) -> None:
    """A properly tightened solver draws no complaint."""
    import logging

    toy_model.tolerance = 1e-9
    with caplog.at_level(logging.WARNING, logger="spectra_cobra.extraction"):
        spectra_me(toy_model, ["R3"], tol=1e-4, seed=0)

    assert not any("model.tolerance" in r.message for r in caplog.records)


def test_lp_weights_follow_the_reaction_not_the_iteration_order(
    toy_model: Model,
) -> None:
    """A seeded run must weight the same reaction the same way every time.

    The random objective coefficients are drawn as one block and zipped
    against the auxiliary variables positionally, and those variables follow
    the order of the identifiers handed in. Callers pass a set, whose
    iteration order varies with ``PYTHONHASHSEED``, so without an explicit
    order the same seed would weight a different reaction in every process
    and the extraction would not be reproducible between runs.
    """
    import numpy as np

    from spectra_cobra._lp import forward_cc
    from spectra_cobra._orientation import reaction_signs

    def coefficients(order: List[str]) -> dict:
        with toy_model:
            forward_cc(
                toy_model,
                order,
                reaction_signs(toy_model),
                1e-3,
                np.random.default_rng(0),
            )
            aux = [v for v in toy_model.variables if v.name.startswith("spectra_aux_")]
            linear = toy_model.objective.get_linear_coefficients(aux)
            return {var.name: float(value) for var, value in linear.items()}

    forwards = coefficients(["R2", "R3", "R4"])
    backwards = coefficients(["R4", "R3", "R2"])
    assert forwards, "the LP should have built auxiliary variables"
    assert forwards == backwards


# --------------------------------------------------------------------------
# What the extracted model carries over.
#
# ``_extract`` builds the result from the kept reactions instead of copying
# the whole model and deleting the rest, because deleting is quadratic in the
# number of removals. Everything the old path preserved by virtue of copying
# now has to be carried deliberately, so each of those things is pinned here.
# --------------------------------------------------------------------------

ALL_PROBLEM_TYPES = (MIN_NET_LP, MIN_NET_MILP, GROWTH_OPTIM, TRADE_OFF)


def _decorate(model: Model) -> Model:
    """Hang metadata on a model so extraction has something to lose."""
    from cobra.core.group import Group

    model.name = "decorated toy"
    for rxn in model.reactions:
        rxn.annotation = {"sbo": "SBO:0000176", "rxn": rxn.id}
        rxn.notes = {"origin": rxn.id}
        rxn.subsystem = "toy subsystem"
    for met in model.metabolites:
        met.annotation = {"chebi": f"CHEBI:{len(met.id)}", "met": met.id}
        met.notes = {"origin": met.id}
    # The toy model carries no genes of its own, so give each reaction a rule
    # before annotating what that creates.
    for index, rxn in enumerate(model.reactions):
        rxn.gene_reaction_rule = f"g{index}"
    for gene in model.genes:
        gene.name = f"{gene.id} product"
        gene.annotation = {"ncbigene": f"gene:{gene.id}", "sbo": "SBO:0000243"}
    model.add_groups(
        [Group("g_all", name="everything", members=list(model.reactions)),
         Group("g_r3", name="just R3", members=[model.reactions.R3]),
         # Every member of this one is dropped when R3 is the core, so the
         # group ends up empty -- which is what deleting the reactions out
         # from under it used to leave behind.
         Group("g_gone", name="dropped branch",
               members=[model.reactions.R7, model.reactions.R8])]
    )
    return model


@pytest.mark.parametrize("problem_type", ALL_PROBLEM_TYPES)
def test_extraction_keeps_annotations_and_notes(
    toy_model: Model, problem_type: str
) -> None:
    """Reaction and metabolite metadata survives every formulation."""
    _decorate(toy_model)
    extracted = spectra_me(
        toy_model, ["R3"], problem_type=problem_type, seed=0
    )
    for rxn in extracted.reactions:
        original = toy_model.reactions.get_by_id(rxn.id)
        assert rxn.annotation == original.annotation, rxn.id
        assert rxn.notes == original.notes, rxn.id
        assert rxn.subsystem == original.subsystem, rxn.id
    for met in extracted.metabolites:
        original = toy_model.metabolites.get_by_id(met.id)
        assert met.annotation == original.annotation, met.id
        assert met.notes == original.notes, met.id


def test_extraction_keeps_groups(toy_model: Model) -> None:
    """Groups come across, with membership cut down to what survived.

    A group whose members were all dropped stays behind as an empty group,
    which is what deleting the reactions out from under it used to do.
    """
    _decorate(toy_model)
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    assert {g.id for g in extracted.groups} == {
        "g_all", "g_r3", "g_gone"}
    assert not extracted.groups.get_by_id("g_gone").members, (
        "a group whose members all went should survive, empty")
    kept = _ids(extracted)
    for group in extracted.groups:
        original = toy_model.groups.get_by_id(group.id)
        expected = {m.id for m in original.members} & kept
        assert {m.id for m in group.members} == expected


def test_extraction_keeps_the_objective(toy_model: Model) -> None:
    """Objective coefficients and direction follow the model."""
    toy_model.objective = {toy_model.reactions.R3: 2.0}
    toy_model.objective.direction = "min"
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    assert extracted.objective.direction == "min"
    assert {
        r.id: r.objective_coefficient
        for r in extracted.reactions
        if r.objective_coefficient
    } == {"R3": 2.0}


def test_objective_direction_survives_losing_every_objective_reaction(
    toy_model: Model,
) -> None:
    """Direction is not collateral damage when the objective empties out.

    Setting it only when some coefficient survived silently reverted a
    minimisation to a maximisation.
    """
    toy_model.objective = {toy_model.reactions.R7: 1.0}
    toy_model.objective.direction = "min"
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    assert "R7" not in _ids(extracted)
    assert extracted.objective.direction == "min"


def test_extraction_keeps_the_tolerance(toy_model: Model) -> None:
    """``model.tolerance`` is not reset to the cobrapy default."""
    toy_model.tolerance = 1e-9
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    assert extracted.tolerance == 1e-9


def test_extraction_keeps_source_ordering(toy_model: Model) -> None:
    """Elements keep the order they had, not the order they were built in.

    Content-wise comparisons are blind to this: they key on identifiers. Only
    sequence catches it, and index-based access depends on it.
    """
    # R11 drags in R9 and R10, whose identifiers sort before R5's -- so
    # alphabetical order differs from the model's and the assertion can
    # actually fail. A core of R3 keeps R1..R5, where the two coincide.
    extracted = spectra_me(toy_model, ["R11"], problem_type=MIN_NET_LP, seed=0)
    kept = _ids(extracted)
    assert sorted(kept) != [r.id for r in toy_model.reactions if r.id in kept], (
        "this core must keep reactions whose sorted order differs from the "
        "model's, or the test cannot detect re-ordering"
    )
    assert [r.id for r in extracted.reactions] == [
        r.id for r in toy_model.reactions if r.id in kept
    ]
    met_ids = {m.id for m in extracted.metabolites}
    assert [m.id for m in extracted.metabolites] == [
        m.id for m in toy_model.metabolites if m.id in met_ids
    ]


def test_extracting_everything_reproduces_the_model(toy_model: Model) -> None:
    """Keeping every reaction gives back an identical model, not merely an
    equivalent one -- same contents in the same order."""
    from cobra.io.dict import model_to_dict

    _decorate(toy_model)
    everything = [r.id for r in toy_model.reactions]
    extracted = spectra_me(
        toy_model, everything, problem_type=MIN_NET_LP, remove_genes=False,
        seed=0,
    )
    assert model_to_dict(extracted) == model_to_dict(toy_model)


def test_the_source_model_is_unharmed(toy_model: Model) -> None:
    """Extraction detaches objects from the source to copy them; it has to
    put them back."""
    before = ([r.id for r in toy_model.reactions],
              [m.id for m in toy_model.metabolites])
    spectra_me(toy_model, ["R3"], seed=0)
    assert ([r.id for r in toy_model.reactions],
            [m.id for m in toy_model.metabolites]) == before
    assert all(r.model is toy_model for r in toy_model.reactions)
    assert all(m.model is toy_model for m in toy_model.metabolites)
    assert all(g.model is toy_model for g in toy_model.genes)


def test_ccme_several_solutions_carry_metadata(toy_model: Model) -> None:
    """``spectra_ccme`` extracts its first solution directly and the rest
    through ``spectra_me``; both paths have to preserve metadata.

    The toy model is used rather than the blocked one because a second
    solution has to exist: its three parallel routes to D give the exclusion
    something to find.
    """
    _decorate(toy_model)
    models, _dropped = spectra_ccme(
        toy_model, ["R5"], n_solutions=2,
        alt_solution_method=PATHWAY_EXCLUSION, problem_type=MIN_NET_MILP,
        seed=0,
    )
    assert isinstance(models, list) and len(models) == 2
    for extracted in models:
        assert extracted.name == "decorated toy"
        for rxn in extracted.reactions:
            assert rxn.annotation == (
                toy_model.reactions.get_by_id(rxn.id).annotation
            )


def test_extraction_keeps_gene_annotations(toy_model: Model) -> None:
    """Gene metadata survives, not merely the gene identifiers.

    ``Reaction.__getstate__`` serialises the GPR to a string, so the rule is
    re-parsed when the reactions are added and the genes arrive as fresh
    objects holding nothing but an identifier. Without putting their metadata
    back, a model whose genes are annotated loses every annotation: iJO1366
    has 1366 of them. Recon3D and Harvey have none, which is why comparing
    extractions of those models said nothing about it.
    """
    _decorate(toy_model)
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    assert extracted.genes, "the toy model should carry genes"
    for gene in extracted.genes:
        original = toy_model.genes.get_by_id(gene.id)
        assert gene.annotation == original.annotation, gene.id
        assert gene.name == original.name, gene.id


def test_extraction_keeps_gene_annotations_when_genes_are_kept(
    toy_model: Model,
) -> None:
    """The same, for the genes put back by ``remove_genes=False``."""
    _decorate(toy_model)
    extracted = spectra_me(toy_model, ["R3"], remove_genes=False, seed=0)
    assert {g.id for g in extracted.genes} == {g.id for g in toy_model.genes}
    for gene in extracted.genes:
        original = toy_model.genes.get_by_id(gene.id)
        assert gene.annotation == original.annotation, gene.id
