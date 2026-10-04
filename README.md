# SPECTRA for cobrapy

**S**calable **P**latform for **E**xtracting **C**onstraint-based **T**op-down
**R**econstructions and **A**nalysis, built on
[cobrapy](https://github.com/opencobra/cobrapy).

SPECTRA reconstructs metabolic networks from multi-omics data at a range of
biological scales: minimal reactomes, context-specific models, gap-filled
reconstructions, minimal microbiomes, microbial community models and
multi-tissue models. What changes between them is the universal model, the
evidence supplied and the objective chosen, not the routine called.

SPECTRA is also available as a
[MATLAB package](https://github.com/NiravBhattLab/SPECTRA) for the COBRA
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
| `"tradeOff"`    | maximise weighted count of kept reactions| **any real**  | 1 MILP |
| `"growthOptim"` | maximise biomass minus weighted flux     | non-negative  | 1 LP   |

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

## Communities and minimal microbiomes

`build_community_model` joins organism models into a community trading
through a shared pool, with every reaction coupled to its organism's biomass
so an absent organism carries no flux at all.

```python
from spectra_cobra import build_community_model, minimal_microbiome

community = build_community_model([a, b, c], organisms=["A", "B", "C"])
result = minimal_microbiome(community, products=["EX_but_u"])
print(result.present, result.membership)
```

`minimal_microbiome` is the `minNetMILP` formulation pointed at organisms:
each biomass reaction gets a binary through `indicator_reactions` and a
weight of 1, everything else is weighted 0, so minimising the weighted count
minimises the membership vector. Growth and production requirements are
measured on the full community first, then imposed as a fraction of it.

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

## Relation to cobrapy's `spectra_cc`

The consistency check is also being contributed to cobrapy itself as
`cobra.flux_analysis.spectra_cc`. The two differ deliberately:

|                     | cobrapy's            | this package's          |
|---------------------|----------------------|-------------------------|
| `consistency_type`  | steady state only    | both modes              |
| detection default   | `model.tolerance`    | `0.99 * tol`            |

cobrapy's version counts any nonzero flux, so it agrees exactly with `fastcc`
and `find_blocked_reactions`. This one is stricter by default: a reaction must
reach `tol` itself — which also drops reactions that can carry *some* flux but
never as much as `tol`. Either behaviour is reachable from either package by
setting the cutoff explicitly.

## Validation

### Parity with the published results

[`tests/test_toy_model_parity.py`](tests/test_toy_model_parity.py) rebuilds the
published toy models and experiments and checks this implementation lands on
the same answers. Expected values are derived from each model's own
stoichiometry rather than copied from a published run, so they are checked
rather than assumed.

**Consistency, topology versus stoichiometry** — consistent reactions found:

| model | stoichiometry | topology |
|---|---|---|
| `topology_toy_model(n=1)` | 0 of 6 | 6 of 6 |
| `topology_toy_model(n=2)` | 6 of 6 | 6 of 6 |
| `topology_toy_model(n=3)` | 0 of 6 | 6 of 6 |
| `get_cc_toy_model_1` | 7 of 8 (`r3` blocked) | 8 of 8 |
| `get_cc_toy_model_2` | 7 of 8 (`r3` blocked) | 7 of 8 (`r3` blocked) |

The `topology_toy_model` rows are the result the topology mode exists for: at
steady state the mass balances force `v2 · (n − 2) = 0`, so only `n = 2`
carries any flux at all, while the accumulation condition admits every `n`.
The two `get_cc_toy_model` rows are the instructive contrast — model 1's `C` is
produced and never consumed, which accumulation rescues; model 2's `C` is
consumed and never produced, which it cannot, since `S·v ≥ 0` lets a
metabolite pile up but not appear from nothing.

**Extraction, topology versus stoichiometry** — with the T1 export
as the sole core reaction reproduces the same split: under stoichiometry the
core is blocked and reported for `n = 1` and `n = 3`, while topology recovers
the full six-reaction network for every `n`.

**Objectives on the three-pathway model** — the formulations on
`three_pathway_toy_model`, with the published weights:

| `problem_type` | result | reactions |
|---|---|---|
| `minNetLP` | `r1 r2 r3 r4 r5` | 5 |
| `minNetMILP` | `r5` + a three-reaction route | 4 |
| `tradeOff` | `r5 r6 r7 r8` | 4 |
| `growthOptim` | everything | 11 |

This model separates the formulations by construction: `r1`–`r4` yields one
`d` per unit of flux over four reactions, while the two shorter routes yield
half a `d` each over three. So the fewest-reaction route is not the
least-flux route, and `minNetLP` (5 reactions) genuinely diverges from
`minNetMILP` (4). `tradeOff` picks `r6`–`r8` because that is the only route
whose published weights sum positive (+2, against −1 each for the others).

Pathway exclusion asked for five solutions returns **exactly the three routes
the network has**, with no duplicates, then stops.

### Genome scale: Recon3D

Note that the published models are *already* consistency checked — `UpdatedRecon3D.mat` and `consRecon3DGeneSymbol.mat` (11303
reactions) and the PCOS study's `Reconmodel.mat` (10600) all have zero blocked
reactions. Running the check on those confirms only that the files are what
they claim. To see it do work, Recon3D was restricted to a defined medium,
which is where context-specific modelling actually starts:

| | reactions | blocked | LPs | time |
|---|---|---|---|---|
| Recon3D, defined medium | 11303 → **8406** | 2897 | 8 | 30 s |
| iJO1366, as distributed | 2583 → **1705** | 878 | 12 | 5 s |

Both derived models are **fixed points** — re-running the check with a
different seed removes nothing — and cobrapy's independent FVA-based
`find_blocked_reactions` reports zero blocked in each.

Where `spectra_cc` and `fastcc` disagree on the defined-medium Recon3D (six
reactions), FVA says `spectra_cc` is right on all six: `fastcc` produced three
false positives and three false negatives. `fastcc` is also not reproducible
at this scale — across five runs it returned 8403, 8405, 8406, 8408 and 8410
on the same model, because its singleton phase draws an arbitrary element from
a Python set, whereas seeded `spectra_cc` returned 8406 every time.

Extraction was then checked with random core sets — twelve trials across four
core sizes, verifying both promises: every core reaction present, and the
result flux consistent. Consistency was checked twice, with this package's
own check and with cobrapy's independent FVA-based `find_blocked_reactions`,
which agreed on every count.

| core | trials | extracted | core present | blocked | of which core |
|---|---|---|---|---|---|
| 10 | 3 | 167–177 | all | 0, 0, 0 | — |
| 50 | 3 | 525–927 | all | 0, 0, **112** | 4 |
| 200 | 3 | 1495–1791 | all (1 failed) | **37**, **19** | 2, 2 |
| 500 | 3 | 2752–2792 | all | **15**, **16**, 0 | 5, 3 |

Core reactions are always present. Whether they can actually *carry flux*
turned out to hinge on two numerical settings rather than on the algorithm:

| Configuration | dead core reactions | trials affected |
|---|---|---|
| inclusion cutoff floored at `model.tolerance`, solver `1e-7` | 16 | 5 of 12, + 1 failure |
| cutoff `tol * 1e-7`, solver `1e-7` | 4 | 2 of 12 |
| cutoff `tol * 1e-7`, **solver `1e-9`** | **0** | **none** |

So: do not floor the inclusion cutoff, and **set `model.tolerance = 1e-9`**
(the lowest Gurobi and CPLEX accept) before extracting. The extraction reads
its answer off an LP whose mass balance holds only to the solver's tolerance,
and at the `1e-7` default that residual is enough to break chains carrying
flux of order `tol`. `spectra_me` warns when the two are within a factor of
1e4.

The cause is the inclusion rule: keep a reaction if its
extraction-LP flux exceeds `tol * 1e-7`, floored here at `model.tolerance`, so
around `1e-7`. A handful of reactions end up carrying flux *just under* that
cutoff while doing load-bearing balancing work. Dropping them leaves a
mass-balance residual of the same order — in the trial examined, up to
`3.5e-07` on `coa[c]`, on 3 of 972 metabolites — which is above the solver's
feasibility tolerance, so the kept flux vector is not actually feasible in the
extracted model. For the chains that relied on those sub-cutoff reactions the
only feasible flux is then exactly zero: all 25 sampled blocked reactions have
an attainable flux of **0** in the extracted model, despite carrying ~`1e-4`
in the LP. They are mostly exchange and transport pairs.

A blocked *core* reaction is the damaging case: it is present, so the model
looks right, but it cannot play the role it was chosen for. Use
`check_extraction` to catch it:

```python
from spectra_cobra import check_extraction

report = check_extraction(extracted, core, tol=1e-4)
assert report.is_valid, report.summary()
```

`minNetMILP` should not have this failure mode: a binary at zero forces its
reaction's flux to *exactly* zero, so nothing sub-cutoff is discarded and the
kept flux vector is exactly balanced. That is an argument from the
formulation, not a measurement — the genome-scale MILP did not finish within
15 minutes on Gurobi at this core size, so it is untested here and the cost is
real. Pruning the blocked reactions afterwards does yield a consistent model,
but it removed two to five core reactions in these trials.

One of the twelve trials raised `SpectraError` because the direction phase's
convex combination cancelled out on two core reactions; a different `seed`
fixes it.

### Other checks

The test suite is 63 tests run against each available solver. Beyond that,
the port has been checked end-to-end on the models bundled with cobrapy,
with glpk and Gurobi agreeing on every result:

| | textbook (95 rxns) | iJO1366 (2583 rxns) |
|---|---|---|
| `spectra_cc` | 87 consistent | 1705 consistent, 10 LPs, 4.9 s |
| `minNetLP` | 27 rxns | 81 rxns, 2.4 s |
| `minNetMILP` | 16 rxns | 61 rxns, 69 s |
| `growthOptim` | 27 rxns | 83 rxns, 2.4 s |
| `tradeOff` | 87 rxns | 1549 rxns, 3.4 s |

Two things worth reading off that table. The consistency check agrees exactly
with cobrapy's own `fastcc` and `find_blocked_reactions` on both models. And
`minNetMILP` lands below `minNetLP` on both, as it should: the former
minimises the reaction count exactly where the latter minimises its L1
relaxation. `tradeOff` keeps nearly everything here only because the default
weights are all +1, which rewards including every reaction; it is meant to be
used with weights that carry both signs.

## Citation

> S, P. K., Sridhar, S., Alsmadi, N., Mahadevan, R., & Bhatt, N. P. (2026).
> *Generalist method to reconstruct metabolic networks from multi-omics data
> at large-scale.* bioRxiv. https://doi.org/10.64898/2026.04.02.716249

## License

MIT. See [LICENSE](LICENSE).
