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
    kept = min_net_milp(
        toy_model, _all_free(toy_model), _ones(toy_model), 1e-4
    ).included
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
    assert (
        min_net_lp(toy_model, directions, weights, 1e-4)
        == min_net_milp(toy_model, directions, weights, 1e-4).included
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

    kept = trade_off(toy_model, _all_free(toy_model), weights, 1e-4).included
    assert {"R1", "R2", "R3", "R4", "R5"} <= kept
    assert not {"R7", "R10"} & kept


def test_trade_off_with_all_negative_weights_keeps_nothing(
    toy_model: Model,
) -> None:
    """If every reaction is penalised, the empty network wins."""
    weights = {r.id: -1.0 for r in toy_model.reactions}
    assert trade_off(toy_model, _all_free(toy_model), weights, 1e-4).included == set()


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


def test_inclusion_cutoff_defaults_below_solver_tolerance(toy_model: Model) -> None:
    """The default cutoff is tol * 1e-7, not floored at model.tolerance.

    Flooring it at the solver tolerance is what left extracted models holding
    blocked core reactions, so the floor's absence is the behaviour under
    test rather than an oversight.
    """
    import spectra_cobra.formulations as module
    from spectra_cobra.formulations import INCLUSION_CUTOFF_FACTOR

    tol = 1e-4
    toy_model.tolerance = 1e-7
    assert tol * INCLUSION_CUTOFF_FACTOR < toy_model.tolerance, "floor would bite"

    directions = _all_free(toy_model)
    directions["R3"] = 1

    captured = {}
    original = module._included_reactions

    def spy(model, t, inclusion_cutoff=None):
        captured["cutoff"] = (
            inclusion_cutoff
            if inclusion_cutoff is not None
            else t * INCLUSION_CUTOFF_FACTOR
        )
        return original(model, t, inclusion_cutoff)

    module._included_reactions = spy
    try:
        min_net_lp(toy_model, directions, _ones(toy_model), tol)
    finally:
        module._included_reactions = original

    assert captured["cutoff"] == tol * INCLUSION_CUTOFF_FACTOR


def test_raising_the_inclusion_cutoff_shrinks_the_model(toy_model: Model) -> None:
    """A cutoff above every flux keeps only what the directions force."""
    directions = _all_free(toy_model)
    directions["R3"] = 1

    default = min_net_lp(toy_model, directions, _ones(toy_model), 1e-4)
    strict = min_net_lp(
        toy_model, directions, _ones(toy_model), 1e-4, inclusion_cutoff=1e3
    )

    assert strict <= default
    assert len(strict) < len(default)


def test_inclusion_cutoff_reaches_spectra_me(toy_model: Model) -> None:
    """spectra_me passes the cutoff down to the LP formulations."""
    from spectra_cobra import spectra_me

    default = spectra_me(toy_model, ["R3"], tol=1e-4, seed=0)
    strict = spectra_me(toy_model, ["R3"], tol=1e-4, seed=0, inclusion_cutoff=1e3)

    assert {r.id for r in strict.reactions} < {r.id for r in default.reactions}


def test_a_milp_keeps_what_carries_flux_whatever_its_binary_says(
    toy_model: Model,
) -> None:
    """A binary at zero still permits a trace flux, so flux decides too.

    ``v_i <= ub_i z_i`` with a bound of 1000 leaves room for ``|v_i|`` up to
    ``1000 * IntFeasTol`` even at ``z_i = 0``. On iJO1366 that was enough to
    drop nineteen reactions carrying 2e-07 apiece, which the solution needed,
    leaving an extracted model that could not grow at all. So the kept set is
    the union of the selected binaries and whatever carries flux, while
    ``selected`` still reports the binaries for the exclusion constraints.
    """
    directions = _all_free(toy_model)
    directions["R3"] = 1
    weights = _ones(toy_model)

    solution = min_net_milp(toy_model, directions, weights, 1e-4)

    assert solution.selected <= solution.included
    assert "R3" in solution.included


def test_milp_solution_unions_the_flux_with_the_binaries(toy_model: Model) -> None:
    """With no binaries at all, everything kept has to come from the flux."""
    from spectra_cobra.formulations import _milp_solution

    with toy_model:
        toy_model.objective = toy_model.reactions.R5
        toy_model.slim_optimize()
        solution = _milp_solution(toy_model, {}, {}, 1e-4)
        carrying = {rxn.id for rxn in toy_model.reactions if abs(rxn.flux) > 1e-9}

    assert solution.selected == set(), "no indicators means nothing selected"
    assert carrying, "the toy model should be carrying some flux here"
    assert carrying <= solution.included


class _StubSolver:
    """A solver that reports a status and either holds a solution or does not."""

    class _Objective:
        def __init__(self, value):
            self._value = value

        @property
        def value(self):
            if isinstance(self._value, Exception):
                raise self._value
            return self._value

    def __init__(self, status, value):
        self.status = status
        self.objective = self._Objective(value)


class _StubModel:
    """Just enough of a model for :func:`_check_status` to inspect."""

    def __init__(self, status, value):
        self.solver = _StubSolver(status, value)


def test_a_time_limited_milp_keeps_its_incumbent(caplog) -> None:
    """Stopping at the time limit with a solution is the point of the limit.

    Gurobi reports ``time_limit`` rather than ``feasible`` when it stops
    holding an incumbent, so leaving that status out made ``time_limit``
    useless: the solution it had paid for was discarded.
    """
    from spectra_cobra.formulations import _check_status

    model = _StubModel("time_limit", 90.0)
    with caplog.at_level("WARNING"):
        assert _check_status(model, "minNetMILP problem", milp=True) == "time_limit"
    assert "rather than optimality" in caplog.text


def test_a_time_limited_milp_without_a_solution_is_an_error() -> None:
    """The same status with nothing to read has to be reported, not returned."""
    from spectra_cobra.exceptions import SpectraSolverError
    from spectra_cobra.formulations import _check_status

    for value in (None, float("nan"), RuntimeError("no solution available")):
        model = _StubModel("time_limit", value)
        with pytest.raises(SpectraSolverError, match="before finding any feasible"):
            _check_status(model, "minNetMILP problem", milp=True)


def test_an_lp_still_demands_optimality() -> None:
    """The relaxed statuses are for mixed-integer solves only."""
    from spectra_cobra.exceptions import SpectraSolverError
    from spectra_cobra.formulations import _check_status

    model = _StubModel("time_limit", 90.0)
    with pytest.raises(SpectraSolverError, match="instead of reaching optimality"):
        _check_status(model, "minNetLP problem")


def test_indicator_reactions_restricts_the_binaries(toy_model: Model) -> None:
    """Only the named reactions get a binary; the rest stay continuous."""
    directions = _all_free(toy_model)
    directions["R3"] = 1
    seen = []

    original = Model.slim_optimize

    def spy(self, *args, **kwargs):
        seen.append({v.name for v in self.variables if v.name.startswith("spectra_z_")})
        return original(self, *args, **kwargs)

    Model.slim_optimize = spy
    try:
        min_net_milp(
            toy_model,
            directions,
            _ones(toy_model),
            1e-4,
            indicator_reactions=["R7", "R10"],
        )
    finally:
        Model.slim_optimize = original

    assert seen[0] == {"spectra_z_R7", "spectra_z_R10"}


def test_reactions_without_an_indicator_are_always_kept(toy_model: Model) -> None:
    """Nothing is choosing over them, so nothing should discard them."""
    directions = _all_free(toy_model)
    directions["R3"] = 1

    solution = min_net_milp(
        toy_model,
        directions,
        _ones(toy_model),
        1e-4,
        indicator_reactions=["R7"],
    )

    # Everything free except R7 is outside the selection and so retained.
    free_but_unselected = {
        rxn.id for rxn in toy_model.reactions if rxn.id not in ("R3", "R7")
    }
    assert free_but_unselected <= solution.included
    assert "R3" in solution.included


def test_no_indicator_reactions_is_the_old_behaviour(toy_model: Model) -> None:
    """Leaving it unset must select over every free reaction as before."""
    directions = _all_free(toy_model)
    directions["R3"] = 1
    weights = _ones(toy_model)

    everything = min_net_milp(toy_model, directions, weights, 1e-4)
    explicit = min_net_milp(
        toy_model,
        directions,
        weights,
        1e-4,
        indicator_reactions=[r.id for r in toy_model.reactions if r.id != "R3"],
    )

    assert everything.included == explicit.included


def test_indicator_reactions_rejects_a_directed_reaction(toy_model: Model) -> None:
    """A reaction already forced in cannot also be selected over."""
    directions = _all_free(toy_model)
    directions["R3"] = 1

    with pytest.raises(SpectraError, match="must be free reactions"):
        min_net_milp(
            toy_model,
            directions,
            _ones(toy_model),
            1e-4,
            indicator_reactions=["R3"],
        )


def test_trade_off_also_takes_indicator_reactions(toy_model: Model) -> None:
    """The option is on both mixed-integer formulations, not just one."""
    weights = {r.id: -1.0 for r in toy_model.reactions}
    weights["R7"] = 10.0

    solution = trade_off(
        toy_model, _all_free(toy_model), weights, 1e-4, indicator_reactions=["R7"]
    )

    assert solution.selected <= {"R7"}
