"""Shared fixtures for the SPECTRA test suite."""

from typing import List

import pytest
from cobra import Model, Reaction
from cobra.util.solver import solvers

#: Solvers that cannot run the whole suite, so are left out of the matrix.
#: ``glpk_exact`` is an exact-arithmetic LP solver and ``scipy`` a thin LP
#: wrapper; neither supports the mixed-integer problems the extraction
#: formulations need.
UNSUPPORTED_SOLVERS = {"glpk_exact", "scipy"}

#: Every solver cobrapy can see that can run the suite, so that each test runs
#: on all of them.
AVAILABLE_SOLVERS = sorted(set(solvers) - UNSUPPORTED_SOLVERS)


@pytest.fixture(params=AVAILABLE_SOLVERS, scope="session")
def solver(request) -> str:
    """Return each available solver in turn.

    Every model fixture sets this on the model it builds, so the whole suite
    runs once per solver rather than only on whichever one happens to be the
    default.

    Parameters
    ----------
    request : pytest.FixtureRequest
        The fixture request carrying the solver name.

    Returns
    -------
    str
        The solver's name.

    """
    return request.param


def _build(name: str, reactions: List[tuple], solver: str) -> Model:
    """Build a toy model from a list of reaction specifications.

    Parameters
    ----------
    name : str
        The model's name.
    reactions : list of tuple
        One ``(id, reaction_string, lower_bound, upper_bound)`` per reaction.
    solver : str
        The solver to attach to the model.

    Returns
    -------
    cobra.Model
        The assembled model.

    """
    model = Model(name)
    built = []
    for rxn_id, _, _, _ in reactions:
        built.append(Reaction(rxn_id))
    model.add_reactions(built)
    for rxn_id, equation, lower, upper in reactions:
        reaction = model.reactions.get_by_id(rxn_id)
        reaction.reaction = equation
        reaction.bounds = (lower, upper)
    model.solver = solver
    return model


@pytest.fixture(scope="function")
def toy_model(solver: str) -> Model:
    """Return the three-pathway toy model from the SPECTRA README.

    Eight metabolites and eleven reactions, with the media reactions R1, R6
    and R9 constrained to an uptake of 1 and everything else to 10. This is
    the model the README's quick-start walks through.

    Returns
    -------
    cobra.Model
        The toy model.

    """
    model = _build(
        "three pathway toy model",
        [
            ("R1", "--> A", 0.0, 1.0),
            ("R2", "A --> B", 0.0, 10.0),
            ("R3", "B --> C", 0.0, 10.0),
            ("R4", "C --> D", 0.0, 10.0),
            ("R5", "D -->", 0.0, 10.0),
            ("R6", "--> E", 0.0, 1.0),
            ("R7", "E --> F", 0.0, 10.0),
            ("R8", "F --> D", 0.0, 10.0),
            ("R9", "--> G", 0.0, 1.0),
            ("R10", "G --> H", 0.0, 10.0),
            ("R11", "H --> D", 0.0, 10.0),
        ],
        solver,
    )
    model.objective = model.reactions.R5
    return model


@pytest.fixture(scope="function")
def blocked_model(solver: str) -> Model:
    """Return a toy model with a dead end.

    ``R2`` and ``R3`` form a path into ``C``, which nothing consumes, so both
    are blocked at steady state but unblocked under the accumulation
    condition.

    Returns
    -------
    cobra.Model
        The toy model.

    """
    return _build(
        "blocked",
        [
            ("R1", "--> A", 0.0, 10.0),
            ("R2", "A --> B", 0.0, 10.0),
            ("R3", "B --> C", 0.0, 10.0),
            ("R4", "A -->", 0.0, 10.0),
        ],
        solver,
    )


@pytest.fixture(scope="function")
def reversible_model(solver: str) -> Model:
    """Return a toy model whose middle reaction is reversible.

    ``R2`` can run either way, and the only steady state that uses it runs it
    in reverse, so it exercises the reverse LP.

    Returns
    -------
    cobra.Model
        The toy model.

    """
    return _build(
        "reversible",
        [
            ("R1", "--> A", 0.0, 10.0),
            ("R2", "B <=> A", -10.0, 10.0),
            ("R3", "B -->", 0.0, 10.0),
        ],
        solver,
    )


@pytest.fixture(scope="function")
def backwards_model(solver: str) -> Model:
    """Return a toy model whose reactions only carry negative flux.

    Every reaction is capped at zero, so read in reverse ``R1`` to ``R3``
    form a source-to-sink path and ``R4`` is a dead end. This is the case the
    reaction orientation exists for.

    Returns
    -------
    cobra.Model
        The toy model.

    """
    return _build(
        "backwards",
        [
            ("R1", "A -->", -10.0, 0.0),
            ("R2", "B --> A", -10.0, 0.0),
            ("R3", "--> B", -10.0, 0.0),
            ("R4", "C -->", -10.0, 0.0),
        ],
        solver,
    )
