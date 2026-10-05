"""Test copying a model whose solver cannot serialise itself.

``St`` begins the constraint section of an LP file. A reaction of that
name makes optlang's Gurobi copy write a file Gurobi then refuses to
parse, and the model comes back as a ``GurobiError`` that says nothing
about identifiers. The CarveMe universal reconstruction has exactly such
a reaction, for sulfur diffusion, so this is not a hypothetical.
"""

import pytest
from cobra import Metabolite, Model, Reaction

from spectra_cobra import spectra_me
from spectra_cobra._copy import copy_model


def _named(rxn_id: str, solver: str) -> Model:
    """Return a two-reaction model whose producing reaction has this name."""
    model = Model("tiny")
    inside = Metabolite("a_c", compartment="c")
    outside = Metabolite("a_e", compartment="e")
    model.add_metabolites([inside, outside])
    make = Reaction(rxn_id, lower_bound=0.0, upper_bound=1000.0)
    out = Reaction("OUT", lower_bound=0.0, upper_bound=1000.0)
    model.add_reactions([make, out])
    make.add_metabolites({outside: -1.0, inside: 1.0})
    out.add_metabolites({inside: -1.0})
    model.add_boundary(outside, type="exchange")
    model.reactions.get_by_id("EX_a_e").lower_bound = -10.0
    model.objective = out
    model.solver = solver
    return model


@pytest.mark.parametrize("rxn_id", ["St", "ordinary"])
def test_a_model_can_be_copied_whatever_its_reactions_are_called(
    solver: str, rxn_id: str
) -> None:
    """An identifier is not the solver's business."""
    model = _named(rxn_id, solver)

    copied = copy_model(model)

    assert {r.id for r in copied.reactions} == {r.id for r in model.reactions}
    assert {m.id for m in copied.metabolites} == {m.id for m in model.metabolites}
    assert all(
        copied.reactions.get_by_id(r.id).bounds == r.bounds for r in model.reactions
    )
    assert copied.slim_optimize() == pytest.approx(model.slim_optimize())
    copied.reactions.OUT.upper_bound = 0.0
    assert model.reactions.OUT.upper_bound == 1000.0, "the copy must be independent"


def test_extraction_survives_a_reaction_named_after_an_lp_keyword(
    solver: str,
) -> None:
    """The failure reached users through spectra_me, not through copy."""
    model = _named("St", solver)

    extracted = spectra_me(model, ["St"], tol=1e-4, problem_type="minNetLP")

    assert "St" in extracted.reactions


def test_the_fallback_carries_what_the_solver_would_have(monkeypatch) -> None:
    """Forced, so that it is covered whichever solver the suite is on."""
    from cobra.io import load_model

    model = load_model("textbook")

    def refuse(self):
        raise RuntimeError("solver will not serialise")

    monkeypatch.setattr(Model, "copy", refuse)
    copied = copy_model(model)

    assert len(copied.reactions) == len(model.reactions)
    assert len(copied.metabolites) == len(model.metabolites)
    assert len(copied.genes) == len(model.genes), "genes decide what can be removed"
    assert all(
        str(copied.reactions.get_by_id(r.id).gene_reaction_rule)
        == str(r.gene_reaction_rule)
        for r in model.reactions
    )
    assert {r.id for r in copied.reactions if r.objective_coefficient} == {
        r.id for r in model.reactions if r.objective_coefficient
    }
    assert copied.slim_optimize() == pytest.approx(model.slim_optimize(), rel=1e-6)
