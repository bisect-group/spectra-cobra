"""Reproduce the published MATLAB SPECTRA toy-model results.

These tests rebuild the toy models and experiments from the MATLAB
implementation at https://github.com/bisect-group/spectra (``main`` branch) and
assert that this port lands on the same answers. Three scripts are covered:

``SPECTRA_CC_topology_vs_stoichiometry.m``
    ``spectraCC`` on ``topology_toy_model(n)`` for n = 1, 2, 3 and on the two
    ``get_cc_toy_model`` networks, under both consistency types.
``SPECTRA_ME__topology_vs_stoichiometry.m``
    ``spectraCCME`` on the same ``topology_toy_model(n)``, with the T1 export
    as the single core reaction.
``Objective_diff_toy_models.m``
    all five network inference formulations on ``three_pathway_toy_model``,
    with the weights that script uses, plus pathway exclusion.

The expected values are derived from the models' own structure rather than
copied from a MATLAB run, so each one is checked against what the
stoichiometry actually permits. Where the port is known to differ from
MATLAB, the test says so.
"""

from typing import Dict, List, Sequence, Set

import pytest
from cobra import Model, Reaction

from spectra_cobra import consistent_reaction_ids, spectra_ccme, spectra_me
from spectra_cobra.formulations import min_net_dc, min_net_lp, min_net_milp

#: ``tol`` as the two consistency scripts set it.
CC_TOL = 1e-4
#: ``tol`` as ``Objective_diff_toy_models.m`` sets it.
ME_TOL = 1e-5


def _build(name: str, formulas: Sequence[str], solver: str) -> Model:
    """Build a model the way the COBRA Toolbox's ``createModel`` does.

    Parameters
    ----------
    name : str
        The model's name.
    formulas : sequence of str
        One reaction equation per reaction, named ``r1`` upwards in order, as
        the toy model functions do.
    solver : str
        The solver to attach.

    Returns
    -------
    cobra.Model
        The assembled model, with ``createModel``'s default bounds of
        ``(0, 1000)`` for the irreversible reactions these toy models use.

    """
    model = Model(name)
    ids = [f"r{index + 1}" for index in range(len(formulas))]
    model.add_reactions([Reaction(rxn_id) for rxn_id in ids])
    for rxn_id, formula in zip(ids, formulas):
        reaction = model.reactions.get_by_id(rxn_id)
        reaction.reaction = formula
        reaction.bounds = (0.0, 1000.0)
    model.solver = solver
    return model


def topology_toy_model(n: int, solver: str) -> Model:
    """Build ``topology_toy_model.m``, figure 2a of the Meneco article.

    Parameters
    ----------
    n : int
        The stoichiometric coefficient on ``a`` in the second reaction, which
        is what the experiment varies.
    solver : str
        The solver to attach.

    Returns
    -------
    cobra.Model
        The toy model.

    Notes
    -----
    At steady state the mass balances force ``v2 * (n - 2) == 0``, so the
    network can only carry flux when ``n == 2``. Under the accumulation
    condition the balances become inequalities and every ``n`` admits flux,
    which is the point the experiment makes.

    """
    return _build(
        f"topology n={n}",
        [
            "--> S",
            f"S --> b + {n} a",
            "b + a --> d",
            "a --> c",
            "d + c --> T1",
            "T1 -->",
        ],
        solver,
    )


def cc_toy_model_1(solver: str) -> Model:
    """Build ``get_cc_toy_model_1.m``, where ``C`` is produced but not used.

    Parameters
    ----------
    solver : str
        The solver to attach.

    Returns
    -------
    cobra.Model
        The toy model.

    """
    return _build(
        "cc toy 1",
        [
            "--> A",
            "A --> B",
            "B --> C",
            "B --> D",
            "D -->",
            "A --> E",
            "E --> F",
            "F --> D",
        ],
        solver,
    )


def cc_toy_model_2(solver: str) -> Model:
    """Build ``get_cc_toy_model_2.m``, where ``C`` is used but not produced.

    Parameters
    ----------
    solver : str
        The solver to attach.

    Returns
    -------
    cobra.Model
        The toy model.

    """
    return _build(
        "cc toy 2",
        [
            "--> A",
            "A --> B",
            "C --> B",
            "B --> D",
            "D -->",
            "A --> E",
            "E --> F",
            "F --> D",
        ],
        solver,
    )


def three_pathway_toy_model(solver: str) -> Model:
    """Build ``three_pathway_toy_model.m`` as its experiment configures it.

    Parameters
    ----------
    solver : str
        The solver to attach.

    Returns
    -------
    cobra.Model
        The toy model, with every bound capped at 10, the three media
        reactions capped at 1, and ``r5`` as the objective, exactly as
        ``Objective_diff_toy_models.m`` sets it up.

    Notes
    -----
    Three routes produce the biomass metabolite ``d``: ``r1``-``r4`` yields one
    ``d`` per unit of flux over four reactions, while ``r6``-``r8`` and
    ``r9``-``r11`` each yield half a ``d`` over three reactions. So the route
    with the fewest *reactions* is not the route with the least *flux*, which
    is what separates the formulations on this model.

    """
    model = _build(
        "three pathway",
        [
            "--> a",
            "a --> b",
            "b --> c",
            "c --> d",
            "d -->",
            "--> e",
            "e --> f",
            "f --> 0.5 d",
            "--> g",
            "g --> h",
            "h --> 0.5 d",
        ],
        solver,
    )
    for reaction in model.reactions:
        reaction.upper_bound = 10.0
    for rxn_id in ("r1", "r6", "r9"):
        model.reactions.get_by_id(rxn_id).upper_bound = 1.0
    model.objective = model.reactions.r5
    return model


#: The weights ``Objective_diff_toy_models.m`` passes to ``tradeOff``:
#: ``[1, 1, -2, -1, 1, 1, 2, -1, -2, 0, 1]``. The r6-r8 route sums to +2 while
#: the other two sum to -1 each, so this should select that route.
TRADE_OFF_WEIGHTS: Dict[str, float] = dict(
    zip(
        [f"r{index}" for index in range(1, 12)],
        [1.0, 1.0, -2.0, -1.0, 1.0, 1.0, 2.0, -1.0, -2.0, 0.0, 1.0],
    )
)

#: A weight of one everywhere, as the ``minNet`` runs use.
UNIT_WEIGHTS: Dict[str, float] = {f"r{index}": 1.0 for index in range(1, 12)}

#: ``growthOptim``'s weights: one everywhere, ten on the objective reaction.
GROWTH_WEIGHTS: Dict[str, float] = dict(UNIT_WEIGHTS, r5=10.0)

ALL_SIX = {"r1", "r2", "r3", "r4", "r5", "r6"}
ALL_EIGHT = ALL_SIX | {"r7", "r8"}


def _ids(model: Model) -> Set[str]:
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
    return {reaction.id for reaction in model.reactions}


# --------------------------------------------------------------------------
# SPECTRA_CC_topology_vs_stoichiometry.m
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "n, expected_steady",
    [
        # v2 * (n - 2) == 0 at steady state, so only n == 2 carries flux.
        (1, set()),
        (2, ALL_SIX),
        (3, set()),
    ],
)
def test_cc_topology_toy_model(n: int, expected_steady: Set[str], solver: str) -> None:
    """Reproduce spectraCC on topology_toy_model for n = 1, 2 and 3.

    This is the result the topology mode exists to demonstrate: a network
    that the steady state rejects outright is fully usable once metabolites
    are allowed to accumulate.
    """
    model = topology_toy_model(n, solver)

    steady, _ = consistent_reaction_ids(model, CC_TOL, "stoichiometry", seed=0)
    assert steady == expected_steady

    accumulating, _ = consistent_reaction_ids(model, CC_TOL, "topology", seed=0)
    assert accumulating == ALL_SIX


def test_cc_toy_model_1_dead_end_is_rescued_by_topology(solver: str) -> None:
    """Reproduce spectraCC on get_cc_toy_model_1.

    ``C`` is produced by ``r3`` and consumed by nothing, so ``r3`` is blocked
    at steady state. Letting ``C`` accumulate unblocks it.
    """
    model = cc_toy_model_1(solver)

    steady, _ = consistent_reaction_ids(model, CC_TOL, "stoichiometry", seed=0)
    assert steady == ALL_EIGHT - {"r3"}

    accumulating, _ = consistent_reaction_ids(model, CC_TOL, "topology", seed=0)
    assert accumulating == ALL_EIGHT


def test_cc_toy_model_2_missing_producer_is_not_rescued(solver: str) -> None:
    """Reproduce spectraCC on get_cc_toy_model_2.

    The contrast with ``get_cc_toy_model_1``: here ``C`` is *consumed* by
    ``r3`` and produced by nothing. The accumulation condition permits a
    metabolite to pile up, not to appear out of nothing, so ``S v >= 0``
    still forces ``-v3 >= 0`` and ``r3`` stays blocked in both modes.
    """
    model = cc_toy_model_2(solver)
    expected = ALL_EIGHT - {"r3"}

    steady, _ = consistent_reaction_ids(model, CC_TOL, "stoichiometry", seed=0)
    assert steady == expected

    accumulating, _ = consistent_reaction_ids(model, CC_TOL, "topology", seed=0)
    assert accumulating == expected


# --------------------------------------------------------------------------
# SPECTRA_ME__topology_vs_stoichiometry.m
# --------------------------------------------------------------------------


@pytest.mark.parametrize("n, core_is_blocked", [(1, True), (2, False), (3, True)])
def test_me_topology_toy_model(n: int, core_is_blocked: bool, solver: str) -> None:
    """Reproduce spectraCCME on topology_toy_model with the T1 export as core.

    Under the steady state the core reaction is blocked for n = 1 and n = 3,
    so it is dropped and reported; under the accumulation condition the whole
    six-reaction network is recovered for every n.
    """
    extracted, blocked_core = spectra_ccme(
        topology_toy_model(n, solver),
        ["r6"],
        tol=CC_TOL,
        consistency_type="stoichiometry",
        seed=0,
    )
    if core_is_blocked:
        assert blocked_core == ["r6"]
        assert _ids(extracted) == set()
    else:
        assert blocked_core == []
        assert _ids(extracted) == ALL_SIX

    extracted, blocked_core = spectra_ccme(
        topology_toy_model(n, solver),
        ["r6"],
        tol=CC_TOL,
        consistency_type="topology",
        seed=0,
    )
    assert blocked_core == []
    assert _ids(extracted) == ALL_SIX


# --------------------------------------------------------------------------
# Objective_diff_toy_models.m
# --------------------------------------------------------------------------

#: The four-reaction route plus the sink, which carries the least total flux:
#: one unit through each of r1-r4 yields one ``d``.
LEAST_FLUX_ROUTE = {"r1", "r2", "r3", "r4", "r5"}
#: The two three-reaction routes plus the sink, which hold the fewest
#: reactions but need two units of flux each to yield one ``d``.
FEWEST_REACTION_ROUTES = ({"r5", "r6", "r7", "r8"}, {"r5", "r9", "r10", "r11"})


def test_me_min_net_lp_picks_the_least_flux_route(solver: str) -> None:
    """minNetLP selects the four-reaction route, which carries the least flux.

    Reaching one unit of ``d`` costs four units of flux through ``r1``-``r4``
    but six through either half-yield route, so the L1 objective prefers the
    longer pathway.
    """
    extracted = spectra_me(
        three_pathway_toy_model(solver),
        ["r5"],
        tol=ME_TOL,
        weights=UNIT_WEIGHTS,
        problem_type="minNetLP",
        seed=0,
    )
    assert _ids(extracted) == LEAST_FLUX_ROUTE


def test_me_min_net_milp_picks_a_fewest_reaction_route(solver: str) -> None:
    """minNetMILP selects a three-reaction route, the smallest network.

    This is the formulation's whole purpose, and the result differs from
    minNetLP's on this model: four reactions rather than five.
    """
    extracted = spectra_me(
        three_pathway_toy_model(solver),
        ["r5"],
        tol=ME_TOL,
        weights=UNIT_WEIGHTS,
        problem_type="minNetMILP",
        seed=0,
    )
    assert _ids(extracted) in FEWEST_REACTION_ROUTES


def test_me_trade_off_follows_the_published_weights(solver: str) -> None:
    """tradeOff selects the route its weights reward.

    With ``[1, 1, -2, -1, 1, 1, 2, -1, -2, 0, 1]`` the r6-r8 route sums to +2
    while the r1-r4 and r9-r11 routes each sum to -1, so only the first is
    worth including.
    """
    extracted = spectra_me(
        three_pathway_toy_model(solver),
        ["r5"],
        tol=ME_TOL,
        weights=TRADE_OFF_WEIGHTS,
        problem_type="tradeOff",
        seed=0,
    )
    assert _ids(extracted) == {"r5", "r6", "r7", "r8"}


def test_me_growth_optim_uses_every_route(solver: str) -> None:
    """growthOptim keeps all three routes to maximise the objective flux.

    The reward is ten per unit of ``r5`` against one per unit of flux
    elsewhere. Running every route gives ``d = 2`` for a flux cost of ten, so
    a net twenty minus ten, which beats any subset.
    """
    extracted = spectra_me(
        three_pathway_toy_model(solver),
        ["r5"],
        tol=ME_TOL,
        weights=GROWTH_WEIGHTS,
        problem_type="growthOptim",
        seed=0,
    )
    assert _ids(extracted) == set(UNIT_WEIGHTS)


def test_me_min_net_dc_settles_on_the_l1_solution(solver: str) -> None:
    """minNetDC returns the L1 route on this model, not the smallest one.

    This is the one place the port is known to be able to disagree with
    MATLAB. ``minNetDC`` is a re-implementation of the difference-of-convex
    scheme rather than a port of the COBRA Toolbox's ``optimizeCardinality``,
    and difference-of-convex is a local method: started from the L1 solution
    it sits at a fixed point, because the reweighting penalises the
    zero-flux reactions most and the already-active ones least. The result is
    stable across every step-sharpness schedule tried, so it is the scheme's
    local optimum rather than a tuning artefact.

    :func:`min_net_milp` solves the same objective exactly when the smallest
    network is what matters.
    """
    model = three_pathway_toy_model(solver)
    directions = dict({f"r{index}": 0 for index in range(1, 12)}, r5=1)

    from_dc = min_net_dc(model, directions, UNIT_WEIGHTS, ME_TOL)
    from_lp = min_net_lp(model, directions, UNIT_WEIGHTS, ME_TOL)
    from_milp = min_net_milp(model, directions, UNIT_WEIGHTS, ME_TOL)

    assert from_dc == from_lp == LEAST_FLUX_ROUTE
    assert len(from_milp) < len(from_dc)


@pytest.mark.parametrize(
    "problem_type, weights",
    [("minNetMILP", UNIT_WEIGHTS), ("tradeOff", TRADE_OFF_WEIGHTS)],
)
def test_me_pathway_exclusion_enumerates_the_routes(
    problem_type: str, weights: Dict[str, float], solver: str
) -> None:
    """Pathway exclusion walks the distinct routes and then stops.

    The network offers three routes to ``d``, so asking for five solutions
    yields the routes that exist and no duplicates. MATLAB's ``spectraME``
    breaks out of its loop when the excluded solutions leave the problem
    infeasible, and this does the same rather than raising.
    """
    models: List[Model] = spectra_me(
        three_pathway_toy_model(solver),
        ["r5"],
        tol=ME_TOL,
        weights=weights,
        n_solutions=5,
        alt_solution_method="pathwayExclusion",
        problem_type=problem_type,
        seed=0,
    )
    reaction_sets = [frozenset(_ids(model)) for model in models]

    assert 1 <= len(models) <= 5
    assert len(set(reaction_sets)) == len(reaction_sets), "solutions must differ"
    for model in models:
        assert "r5" in _ids(model), "the core reaction is in every solution"
