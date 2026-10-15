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

| `problem_type`  | Objective                                | Weights       | Solves |
|-----------------|------------------------------------------|---------------|--------|
| `"minNetLP"`    | minimise weighted total absolute flux    | non-negative  | 1 LP   |
| `"minNetMILP"`  | minimise weighted reaction count         | non-negative  | 1 MILP |
| `"tradeOff"`    | maximise summed weight of kept reactions | **any real**  | 1 MILP |
| `"growthOptim"` | maximise biomass minus weighted flux     | non-negative  | 1 LP   |

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

## Gap-filling

`gapfill_for_growth` adds the cheapest set of universal reactions that lets a
draft grow in each of a panel of media; `gapfill_for_tasks` does the same for
a task list.

```python
from spectra_cobra import gapfill_for_growth

M9 = {"EX_pi_e": -1000, "EX_h2o_e": -1000, "EX_h_e": -1000,
      "EX_nh4_e": -1000, "EX_o2_e": -1000, "EX_co2_e": -1000}
media = {c: {**M9, f"EX_{c}_e": -10.0} for c in ("glc__D", "ac", "succ")}

result = gapfill_for_growth(draft, universal, media)   # min_growth=0.1
print(result.summary(), result.added)
```

Growth is demanded through a **lower bound on the biomass reaction**, not
through `tol`. Conflating them also demands that rate of every other core
reaction you supply, and makes the answer move with the tolerance; with the
bound, iJO1366 returns the same 15 reactions at every tolerance from `1e-4`
to `1e-7`. `tradeOff` is refused here: it requires every included reaction to
carry at least `tol`, which a growth solution cannot satisfy.

For tasks, the reactions each task cannot do without become the core set.
That is a necessary condition and not a sufficient one — a task with
alternative routes has no essential reactions — so `gapfill_for_tasks`
follows up task by task on whatever still fails.

## Communities, tissues and minimal microbiomes

`build_community_model` joins models into one network whose parts trade
through a **pool**. The parts may be organisms of a community or tissues of
a body — nothing in the construction distinguishes them. Each unit's own
exchange reactions are replaced, so the only way out is through a pool you
nominate, and each unit's reactions are coupled to an **anchor** so a unit
that is not working carries no flux at all.

```python
from spectra_cobra import build_community_model, minimal_microbiome

community = build_community_model([a, b, c], organisms=["A", "B", "C"])
result = minimal_microbiome(community, products=["EX_but_u"])
print(result.present, result.membership)
```

A body needs several pools, and which of them reach the environment has to
be said: a pool that merely links two tissues must stay closed, or those
two could bypass the blood entirely.

```python
from spectra_cobra import build_multi_tissue_model, gapfill_community

body = build_multi_tissue_model(
    gem,
    tissues=["tis1", "tis2", "tis3"],
    pools={"Bl": ["tis1", "tis2"], "tis2_tis3": ["tis2", "tis3"]},
    environment=["Bl"],
    anchor_reactions={t: "DM_atp_c_" for t in ("tis1", "tis2", "tis3")},
    link_bounds=measured_rates,   # per-tissue uptake, from metabolomics
)
print(body.coupling_summary())    # what each tissue got coupled to
```

`gapfill_community` then fills every unit in one solve, so a gap in one can
be closed by another's secretion instead of by a new reaction. That is why
the joint answer is smaller than filling the members one at a time.

`minimal_microbiome` is `spectra_me` with `problem_type="minNetMILP"`
pointed at organisms: each biomass reaction gets a binary through
`indicator_reactions` and a weight of 1, everything else is weighted 0, so
minimising the weighted count minimises the membership vector. Growth and
production requirements are measured on the full community first, then
imposed as a fraction of it.

Every routine in the package reaches the solver through `spectra_cc`,
`spectra_me` or `spectra_ccme` — there is no second implementation of any
formulation. What changes between a context-specific model, a minimal
reactome, a gap-filled draft and a minimal microbiome is the universal
model, the evidence supplied and the objective chosen.

## Minimal reactomes

The opposite of gap-filling: the smallest set of reactions that still meets
a requirement. `minimal_reactome` is `minNetMILP` with no core set at all —
every reaction gets a binary, the requirements enter as bounds, and the
objective counts reactions.

```python
from spectra_cobra import minimal_reactome

result = minimal_reactome(
    model,
    medium={"glucose": glucose, "acetate": acetate},   # must work in both
    growth_fraction=1.0,              # lose no growth at all
    products={"EX_succ_e": 1.0},      # and still secrete succinate
    tasks=task_list,                  # and still perform these
    keep_reactions=["ATPM"],          # and still carry flux through these
)
print(result.summary())
```

The answer is minimal *for an environment*, which is Burgard's original
point: on iJO1366 at full wild-type growth it is 437 reactions on glucose,
431 on acetate, 433 on glycerol and 434 on succinate — 17% of a
2583-reaction model, about 15 s each, every one verified.

A metabolic task defines its own medium, so it cannot share a flux vector
with the growth requirement; neither can a second medium. Both therefore
become separate **conditions**: the network is replicated once per condition
and one binary per reaction governs every copy, so what is minimised is the
size of the union while each condition is satisfied in its own right. On the
textbook model with one growth condition and three tasks that is 51
reactions, against 58 for the union of four separate solves; on iJO1366 with
three media it is 448 against 458.

`preprocess` is on by default and settles what it can exactly before the
MILP starts — reactions that can carry no flux in any condition are deleted,
reactions some condition is infeasible without are forced in without a
binary. On iJO1366 on glucose that is 2057 blocked and 407 essential,
leaving 119 binaries out of 2583, and the objective is 28 either way. It
costs LPs, so it loses on problems that already close quickly and wins on
the ones that do not: at 50% growth and an equal 1200 s budget it returns
425 reactions against 440 without it.

## Citation

> S, P. K., Sridhar, S., Alsmadi, N., Mahadevan, R., & Bhatt, N. P. (2026).
> *Generalist method to reconstruct metabolic networks from multi-omics data
> at large-scale.* bioRxiv. https://doi.org/10.64898/2026.04.02.716249

## License

MIT. See [LICENSE](LICENSE).
