"""Test the SPECTRA consistency check."""

import pytest
from cobra import Model
from cobra.flux_analysis import find_blocked_reactions

from spectra_cobra import (
    SpectraError,
    blocked_reaction_ids,
    consistent_reaction_ids,
    spectra_cc,
)


def test_toy_model_is_fully_consistent(toy_model: Model) -> None:
    """Every reaction of the three-pathway toy model carries flux."""
    consistent = spectra_cc(toy_model, seed=0)
    assert {r.id for r in consistent.reactions} == {r.id for r in toy_model.reactions}


def test_dead_end_is_blocked(blocked_model: Model) -> None:
    """A path into a metabolite nothing consumes is blocked."""
    consistent = spectra_cc(blocked_model, seed=0)
    assert {r.id for r in consistent.reactions} == {"R1", "R4"}


def test_reverse_only_reaction_is_consistent(reversible_model: Model) -> None:
    """A reversible reaction usable only in reverse is not blocked."""
    consistent = spectra_cc(reversible_model, seed=0)
    assert {r.id for r in consistent.reactions} == {"R1", "R2", "R3"}


def test_negative_flux_reactions(backwards_model: Model) -> None:
    """Reactions restricted to negative flux are oriented correctly."""
    consistent = spectra_cc(backwards_model, seed=0)
    assert {r.id for r in consistent.reactions} == {"R1", "R2", "R3"}


def test_agrees_with_find_blocked_reactions(blocked_model: Model) -> None:
    """SPECTRA and cobrapy's FVA-based check agree on what is blocked."""
    blocked_model.tolerance = 1e-7
    spectra_blocked = set(
        blocked_reaction_ids(blocked_model, tol=1e-3, detection_cutoff=1e-7, seed=0)
    )
    assert spectra_blocked == set(find_blocked_reactions(blocked_model))


def test_topology_is_a_superset_of_stoichiometry(blocked_model: Model) -> None:
    """The accumulation condition never finds fewer consistent reactions."""
    steady, _ = consistent_reaction_ids(
        blocked_model, consistency_type="stoichiometry", seed=0
    )
    accumulating, _ = consistent_reaction_ids(
        blocked_model, consistency_type="topology", seed=0
    )
    assert steady < accumulating
    # The dead end into C is fine once C may accumulate.
    assert accumulating == {"R1", "R2", "R3", "R4"}


def test_is_reproducible_given_a_seed(toy_model: Model) -> None:
    """The same seed gives the same answer, and so does a different one."""
    first, lps = consistent_reaction_ids(toy_model, seed=0)
    repeated, repeated_lps = consistent_reaction_ids(toy_model, seed=0)
    other, _ = consistent_reaction_ids(toy_model, seed=11)
    assert first == repeated
    assert lps == repeated_lps
    assert first == other


def test_lp_count_is_reported(toy_model: Model) -> None:
    """The number of LPs solved is a positive multiple of two."""
    _, n_lps = consistent_reaction_ids(toy_model, seed=0)
    assert n_lps > 0 and n_lps % 2 == 0


def test_rejects_unknown_consistency_type(toy_model: Model) -> None:
    """An unrecognised consistency type is an error, not a silent fallback."""
    with pytest.raises(SpectraError, match="consistency_type"):
        spectra_cc(toy_model, consistency_type="nonsense")


def test_does_not_modify_the_input_model(blocked_model: Model) -> None:
    """The input model keeps its reactions, bounds and constraint senses."""
    before_bounds = {r.id: r.bounds for r in blocked_model.reactions}
    before_constraints = {
        m.id: (
            blocked_model.constraints.get(m.id).lb,
            blocked_model.constraints.get(m.id).ub,
        )
        for m in blocked_model.metabolites
    }

    spectra_cc(blocked_model, consistency_type="topology", seed=0)

    assert {r.id: r.bounds for r in blocked_model.reactions} == before_bounds
    assert {
        m.id: (
            blocked_model.constraints.get(m.id).lb,
            blocked_model.constraints.get(m.id).ub,
        )
        for m in blocked_model.metabolites
    } == before_constraints
    assert len(blocked_model.variables) == 2 * len(blocked_model.reactions)
