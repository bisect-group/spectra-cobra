"""Test the network inference formulations directly."""

import pytest
from cobra import Model

from spectra_cobra import (
    SpectraError,
    flux_reducer,
    min_net_lp,
    min_net_milp,
    trade_off,
)
from spectra_cobra._orientation import reaction_signs


def _all_free(model: Model) -> dict:
    """Return a direction of 0 for every reaction, i.e. nothing is core.

    Parameters
    ----------
    model : cobra.Model
        The model to describe.

    Returns
    -------
    dict of {str: int}
        A direction of 0 for each reaction.

    """
    return {rxn.id: 0 for rxn in model.reactions}


def _ones(model: Model) -> dict:
    """Return a weight of one for every reaction.

    Parameters
    ----------
    model : cobra.Model
        The model to describe.

    Returns
    -------
    dict of {str: float}
        A weight of 1.0 for each reaction.

    """
    return {rxn.id: 1.0 for rxn in model.reactions}


def test_min_net_lp_with_no_core_keeps_nothing(toy_model: Model) -> None:
    """With nothing to keep, minimising flux returns the empty network."""
    kept = min_net_lp(toy_model, _all_free(toy_model), _ones(toy_model), 1e-4)
    assert kept == set()


def test_min_net_milp_with_no_core_keeps_nothing(toy_model: Model) -> None:
    """The MILP also settles on the empty network when nothing is forced."""
    kept = min_net_milp(toy_model, _all_free(toy_model), _ones(toy_model), 1e-4)
    assert kept == set()


def test_min_net_lp_honours_a_forced_direction(toy_model: Model) -> None:
    """A direction of 1 forces positive flux and pulls in a pathway."""
    directions = _all_free(toy_model)
    directions["R3"] = 1
    kept = min_net_lp(toy_model, directions, _ones(toy_model), 1e-4)
    assert kept == {"R1", "R2", "R3", "R4", "R5"}


def test_min_net_milp_matches_min_net_lp_here(toy_model: Model) -> None:
    """On this toy model the exact and relaxed formulations agree."""
    directions = _all_free(toy_model)
    directions["R3"] = 1
    weights = _ones(toy_model)
    assert min_net_lp(toy_model, directions, weights, 1e-4) == min_net_milp(
        toy_model, directions, weights, 1e-4
    )


def test_min_net_milp_uses_binary_variables(toy_model: Model) -> None:
    """The MILP really is mixed-integer, not a silently relaxed LP."""
    directions = _all_free(toy_model)
    directions["R3"] = 1
    seen_types = []

    original = Model.slim_optimize

    def spy(self, *args, **kwargs):
        seen_types.append(
            [v.type for v in self.variables if v.name.startswith("spectra_z_")]
        )
        return original(self, *args, **kwargs)

    Model.slim_optimize = spy
    try:
        min_net_milp(toy_model, directions, _ones(toy_model), 1e-4)
    finally:
        Model.slim_optimize = original

    assert seen_types and seen_types[0]
    assert set(seen_types[0]) == {"binary"}


def test_trade_off_splits_reversible_reactions(reversible_model: Model) -> None:
    """A free reversible reaction gets a direction pair of binaries."""
    directions = _all_free(reversible_model)
    seen = []

    original = Model.slim_optimize

    def spy(self, *args, **kwargs):
        seen.append({v.name for v in self.variables if v.name.startswith("spectra_")})
        return original(self, *args, **kwargs)

    Model.slim_optimize = spy
    try:
        trade_off(reversible_model, directions, _ones(reversible_model), 1e-4)
    finally:
        Model.slim_optimize = original

    assert seen
    # R2 is the reversible one, so it is the one that gets a and b.
    assert "spectra_a_R2" in seen[0]
    assert "spectra_b_R2" in seen[0]
    assert "spectra_a_R1" not in seen[0]


def test_trade_off_keeps_positively_weighted_reactions(toy_model: Model) -> None:
    """tradeOff includes a reaction whose weight rewards inclusion."""
    weights = {r.id: -1.0 for r in toy_model.reactions}
    for rxn_id in ("R1", "R2", "R3", "R4", "R5"):
        weights[rxn_id] = 10.0

    kept = trade_off(toy_model, _all_free(toy_model), weights, 1e-4)
    assert {"R1", "R2", "R3", "R4", "R5"} <= kept
    assert not {"R7", "R10"} & kept


def test_trade_off_with_all_negative_weights_keeps_nothing(
    toy_model: Model,
) -> None:
    """If every reaction is penalised, the empty network wins."""
    weights = {r.id: -1.0 for r in toy_model.reactions}
    assert trade_off(toy_model, _all_free(toy_model), weights, 1e-4) == set()


def test_min_net_lp_rejects_a_blocked_forced_reaction(
    blocked_model: Model,
) -> None:
    """Forcing flux through a blocked reaction makes the LP infeasible."""
    directions = _all_free(blocked_model)
    directions["R3"] = 1
    with pytest.raises(SpectraError):
        min_net_lp(blocked_model, directions, _ones(blocked_model), 1e-4)


def test_topology_mode_admits_the_dead_end(blocked_model: Model) -> None:
    """Under the accumulation condition the dead end can carry flux."""
    directions = _all_free(blocked_model)
    directions["R3"] = 1
    kept = min_net_lp(
        blocked_model, directions, _ones(blocked_model), 1e-4, steady_state=False
    )
    assert {"R1", "R2", "R3"} <= kept


def test_negative_direction_forces_reverse_flux(reversible_model: Model) -> None:
    """A direction of -1 forces negative oriented flux."""
    directions = _all_free(reversible_model)
    directions["R2"] = -1
    kept = min_net_lp(reversible_model, directions, _ones(reversible_model), 1e-4)
    assert "R2" in kept


def test_flux_reducer_returns_a_feasible_steady_state(toy_model: Model) -> None:
    """The reduced flux satisfies the mass balance of every metabolite."""
    fluxes = flux_reducer(toy_model)
    assert set(fluxes) == {r.id for r in toy_model.reactions}

    for metabolite in toy_model.metabolites:
        balance = sum(
            reaction.get_coefficient(metabolite.id) * fluxes[reaction.id]
            for reaction in metabolite.reactions
        )
        assert abs(balance) < 1e-6


def test_flux_reducer_is_sparse(toy_model: Model) -> None:
    """Minimising total flux on a model with no demand gives no flux."""
    fluxes = flux_reducer(toy_model)
    assert all(abs(v) < 1e-6 for v in fluxes.values())


def test_flux_reducer_respects_forced_uptake(toy_model: Model) -> None:
    """With an uptake forced on, the reduced flux routes it to a sink."""
    toy_model.reactions.R1.bounds = (1.0, 1.0)
    fluxes = flux_reducer(toy_model)
    assert abs(fluxes["R1"] - 1.0) < 1e-6
    # The A-B-C-D route is the only way out of A.
    assert abs(fluxes["R5"] - 1.0) < 1e-6


def test_orientation_of_negative_only_reactions(backwards_model: Model) -> None:
    """Reactions capped at zero flux are oriented the other way."""
    signs = reaction_signs(backwards_model)
    assert signs == {"R1": -1.0, "R2": -1.0, "R3": -1.0, "R4": -1.0}
