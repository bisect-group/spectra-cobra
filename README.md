# SPECTRA for cobrapy

**S**calable **P**latform for **E**xtracting **C**onstraint-based **T**op-down
**R**econstructions and **A**nalysis, built on
[cobrapy](https://github.com/opencobra/cobrapy).

SPECTRA checks flux consistency and reconstructs metabolic networks from
multi-omics data at a range of biological scales: minimal reactomes,
context-specific models, gap-filled reconstructions, minimal microbiomes,
microbial community models and multi-tissue models. What changes between
them is the universal model, the evidence supplied and the objective
chosen, not the routine called.

If SPECTRA is useful in your work, a star on the repository helps other
people find it.

SPECTRA is also available as a
[MATLAB package](https://github.com/bisect-group/SPECTRA) for the COBRA
Toolbox.

Documentation: <https://spectra-cobra.readthedocs.io>

## Installation

```bash
git clone https://github.com/bisect-group/spectra-cobra.git
cd spectra-cobra
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
from spectra_cobra import spectra_cc, spectra_me

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
from spectra_cobra import spectra_ccme

extracted, blocked_core = spectra_ccme(model, core_reactions=["PGI", "PFK"], tol=1e-3)
if blocked_core:
    print(f"these core reactions are blocked and were dropped: {blocked_core}")
```

## Choosing a formulation

`problem_type` picks what the extraction optimises. What a weight *means*
depends on the choice, so the two columns go together:

| `problem_type`         | Objective                                      | Weights      | Solves |
|------------------------|------------------------------------------------|--------------|--------|
| `"minNetLP"` (default) | minimise weighted total absolute flux          | non-negative | 1 LP   |
| `"minNetMILP"`         | minimise weighted reaction count               | non-negative | 1 MILP |
| `"tradeOff"`           | maximise net weight: positive in, negative out | **any real** | 1 MILP |
| `"growthOptim"`        | maximise biomass minus weighted flux           | non-negative | 1 LP   |

None of this governs the core reactions. Those are constrained to carry at
least `tol` whichever `problem_type` you pick, so they are kept by the
constraints rather than chosen by the objective. What the choice decides is
which *other* reactions come in around them to make that possible.

The two `minNet` variants target the same thing by different routes.
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
it only works with `"minNetMILP"` or `"tradeOff"`. If the network runs out of
alternatives before `n_solutions` is reached, you get the ones that exist
rather than an error, so check `len()` on the result.

## Steady state or accumulation

Every routine takes `consistency_type`:

- `"stoichiometry"` (default) enforces `S·v = 0`, the usual steady state.
- `"topology"` relaxes it to `S·v ≥ 0`, letting metabolites accumulate. This
  is a weaker requirement, so it never finds fewer consistent reactions.

## Reproducibility

The LPs use randomised objective coefficients to break ties between
reactions, so results vary between runs. Pass `seed` to pin them:

```python
extracted = spectra_me(model, core, seed=0)  # same answer every time
```

## Metabolic tasks

A task states something the network must be able to do: take up these
metabolites, produce those, carry flux through this reaction.

The example below asks whether the model can make ATP aerobically. Glucose,
oxygen, phosphate, water and protons may be taken up; carbon dioxide, water
and protons may leave; and the added equation runs ATP hydrolysis at a rate
of at least 1. The task is met only if the network can regenerate the ATP it
is being forced to spend.

```python
from spectra_cobra import MetabolicTask, TaskEquation, check_tasks

task = MetabolicTask(
    "atp_aerobic",
    inputs={"glc__D_e": (0, 1000), "o2_e": (0, 1000), "pi_e": (0, 1000),
            "h2o_e": (0, 1000), "h_e": (0, 1000)},
    outputs={"co2_e": (0, 1000), "h2o_e": (0, 1000), "h_e": (0, 1000)},
    equations=(TaskEquation({"atp_c": -1, "h2o_c": -1,
                             "adp_c": 1, "pi_c": 1, "h_c": 1}, (1.0, 1000.0)),),
)
for result in check_tasks(model, [task]):
    print(result.task.id, result.feasible, result.ok)
```

Bounds default to `(0, 1000)` — allowed, never required — so something needs
a positive lower bound or the task is met by doing nothing. `parse_task_list`
reads the published tab-separated task lists, binding them to a model by
metabolite name and compartment.

## The rest of the package

Gap-filling (`gapfill_for_growth`, `gapfill_for_tasks`), community and
multi-tissue models (`build_community_model`, `build_multi_tissue_model`,
`gapfill_community`), minimal microbiomes (`minimal_microbiome`) and minimal
reactomes (`minimal_reactome`) all reach the solver through `spectra_cc`,
`spectra_me` or `spectra_ccme` — there is no second implementation of any
formulation. What changes between a context-specific model, a minimal
reactome, a gap-filled draft and a minimal microbiome is the universal model,
the evidence supplied and the objective chosen. The
[documentation](https://spectra-cobra.readthedocs.io) covers each in full,
with tutorials that build them up from toy models.

## Citation

> S, P. K., Sridhar, S., Alsmadi, N., Mahadevan, R., & Bhatt, N. P. (2026).
> *Generalist method to reconstruct metabolic networks from multi-omics data
> at large-scale.* bioRxiv. https://doi.org/10.64898/2026.04.02.716249

## License

MIT. See [LICENSE](LICENSE).
