# SPECTRA for cobrapy

**S**calable **P**latform for **E**xtracting **C**onstraint-based **T**op-down
**R**econstructions and **A**nalysis — a Python port of the
[SPECTRA](https://github.com/NiravBhattLab/SPECTRA) MATLAB package, built on
[cobrapy](https://github.com/opencobra/cobrapy) instead of the COBRA Toolbox.

SPECTRA reconstructs metabolic models from multi-omics data, extracting
networks at a range of scales: minimal reactomes, context-specific models,
microbial community models and multi-tissue models.

## Installation

```bash
pip install -e .
```

Needs Python 3.9+ and cobrapy 0.29+. The LP formulations run on any solver
cobrapy supports, including the bundled glpk. The mixed-integer formulations
(`minNetMILP`, `tradeOff`) also work on glpk, but are much faster with Gurobi
or CPLEX:

```bash
pip install -e ".[gurobi]"   # or ".[cplex]"
```

The test suite runs against every solver cobrapy can see, so whichever you
install is covered.

## Quick start

```python
from cobra.io import load_model
from spectra_cobrapy import spectra_cc, spectra_me

model = load_model("textbook")

# 1. Drop the blocked reactions.
consistent = spectra_cc(model, tol=1e-3)

# 2. Extract a model around the reactions you care about.
extracted = spectra_me(consistent, core_reactions=["PGI", "PFK", "FBA"], tol=1e-3)
```

If the universal model is not already flux consistent, skip step 1 and let
`spectra_ccme` do both at once. It reports the core reactions it had to drop
rather than failing:

```python
from spectra_cobrapy import spectra_ccme

extracted, blocked_core = spectra_ccme(model, core_reactions=["PGI", "PFK"], tol=1e-3)
if blocked_core:
    print(f"these core reactions are blocked and were dropped: {blocked_core}")
```

## Choosing a formulation

`problem_type` picks what the extraction optimises. What a weight *means*
depends on the choice, so the two columns go together:

| `problem_type`  | Objective                                | Weights       | Solves |
|-----------------|------------------------------------------|---------------|--------|
| `"minNetLP"`    | minimise weighted total absolute flux    | non-negative  | 1 LP   |
| `"minNetMILP"`  | minimise weighted reaction count         | non-negative  | 1 MILP |
| `"minNetDC"`    | minimise weighted reaction count         | non-negative  | n LPs  |
| `"tradeOff"`    | maximise weighted count of kept reactions| **any real**  | 1 MILP |
| `"growthOptim"` | maximise biomass minus weighted flux     | non-negative  | 1 LP   |

The three `minNet` variants target the same thing by different routes.
`minNetLP` is the L1 relaxation and is by far the cheapest, but it favours
many small fluxes over few large ones, so it tends to keep more reactions
than necessary. `minNetMILP` minimises the count exactly. On the textbook
model the difference is visible: `minNetLP` keeps 27 reactions where
`minNetMILP` keeps 16.

`tradeOff` is the one to use when your omics data gives both positive and
negative evidence, since a positive weight pushes a reaction in and a
negative one pushes it out:

```python
weights = {"R_wanted": 2.0, "R_unwanted": -1.0}  # unlisted reactions default to 1.0
extracted = spectra_me(model, core, weights=weights, problem_type="tradeOff")
```

## Alternative solutions

Ask for several models with `n_solutions`:

```python
# Re-run the direction phase with the LPs in a different order.
models = spectra_me(model, core, n_solutions=5)

# Or forbid each solution already found. Needs a mixed-integer problem_type.
models = spectra_me(
    model, core, n_solutions=5,
    alt_solution_method="pathwayExclusion", problem_type="tradeOff",
)
```

`"coreDirection"` varies which direction the core reactions run in, so the
alternatives may coincide. `"pathwayExclusion"` adds a constraint ruling out
each reaction set already returned, so every model is genuinely distinct, but
it only works with `"minNetMILP"` or `"tradeOff"`.

## Steady state or accumulation

Every routine takes `consistency_type`:

- `"stoichiometry"` (default) enforces `S·v = 0`, the usual steady state.
- `"topology"` relaxes it to `S·v ≥ 0`, letting metabolites accumulate. This
  is a weaker requirement, so it never finds fewer consistent reactions.

## Reproducibility

The LPs use randomised objective coefficients (`unifrnd(1, 1.1)` in MATLAB) to
break ties between reactions, so results vary between runs. Pass `seed` to pin
them:

```python
extracted = spectra_me(model, core, seed=0)  # same answer every time
```

## Relation to cobrapy's `spectra_cc`

The consistency check is also being contributed to cobrapy itself as
`cobra.flux_analysis.spectra_cc`. The two differ deliberately:

|                     | cobrapy's            | this package's          |
|---------------------|----------------------|-------------------------|
| `consistency_type`  | steady state only    | both modes              |
| detection default   | `model.tolerance`    | `0.99 * tol` (MATLAB's) |

cobrapy's version counts any nonzero flux, so it agrees exactly with `fastcc`
and `find_blocked_reactions`. This one defaults to the stricter MATLAB
criterion, where a reaction must reach `tol` itself — which also drops
reactions that can carry *some* flux but never as much as `tol`. Either
behaviour is reachable from either package by setting the cutoff explicitly.

## Differences from the MATLAB implementation

The formulations are ported as-is; the mechanics around them are not.

- **No model mutation.** MATLAB negates the stoichiometric columns of
  reactions with `ub <= 0` and rewrites their bounds. Here that orientation is
  a per-reaction sign applied to the reaction's flux expression, so the input
  model is never modified. The two are equivalent because negating a column
  and negating its variable cancel out.
- **Failures raise.** MATLAB warns and returns `NaN` or an empty solution when
  an LP does not reach optimality. Its callers then misread that: in
  `spectraME`, `abs(NaN) >= tol` is false, so the loop waiting for core
  reactions to be explained never terminates. Here an infeasible core set
  raises `SpectraInfeasibleCoreError` instead of hanging.
- **Two MATLAB bugs fixed.** `spectraCCME.m` passes a misspelled `steadyState`
  to its `growthOptim` and `tradeOff` branches where the variable is
  `steadystate`, so in MATLAB those two problem types raise an
  undefined-variable error rather than running. Both work here.
- **`minNetDC` is a re-implementation, not a port.** MATLAB delegates to the
  COBRA Toolbox's `optimizeCardinality`, which cobrapy has no equivalent of.
  This version minimises the same objective by the standard
  difference-of-convex scheme — approximating the step function with
  `1 - exp(-θ|v|)`, linearising it at the current point, and re-solving while
  increasing `θ`. It targets the same thing but the numbers may differ from
  MATLAB's. Use `minNetMILP` when you need the objective solved exactly.
- **`spectraME2` is not ported.** It is an older LP-only subset of
  `spectraME`, which covers everything it does.

## Validation

The test suite is 48 tests run against each available solver. Beyond that,
the port has been checked end-to-end on the models bundled with cobrapy,
with glpk and Gurobi agreeing on every result:

| | textbook (95 rxns) | iJO1366 (2583 rxns) |
|---|---|---|
| `spectra_cc` | 87 consistent | 1705 consistent, 10 LPs, 4.9 s |
| `minNetLP` | 27 rxns | 81 rxns, 2.4 s |
| `minNetDC` | 27 rxns | 79 rxns, 2.6 s |
| `minNetMILP` | 16 rxns | 61 rxns, 69 s |
| `growthOptim` | 27 rxns | 83 rxns, 2.4 s |
| `tradeOff` | 87 rxns | 1549 rxns, 3.4 s |

Two things worth reading off that table. The consistency check agrees exactly
with cobrapy's own `fastcc` and `find_blocked_reactions` on both models.
And the three `minNet` variants land in the order theory predicts —
`minNetMILP` (exact) below `minNetDC` (cardinality, approximated) below
`minNetLP` (L1 relaxation) — which is the evidence that the DC
re-implementation is doing real work rather than reducing to the L1 problem.
`tradeOff` keeps nearly everything here only because the default weights are
all +1, which rewards including every reaction; it is meant to be used with
weights that carry both signs.

## Citation

> S, P. K., Sridhar, S., Alsmadi, N., Mahadevan, R., & Bhatt, N. P. (2026).
> *Generalist method to reconstruct metabolic networks from multi-omics data
> at large-scale.* bioRxiv. https://doi.org/10.64898/2026.04.02.716249

## License

MIT, as with the MATLAB implementation. See [LICENSE](LICENSE).
