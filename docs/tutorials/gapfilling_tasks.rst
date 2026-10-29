Gap-filling towards metabolic tasks
===================================

A draft reconstruction usually has holes: pathways the organism certainly
runs that the automated reconstruction never assembled. Gap-filling closes
them by borrowing reactions from a universal model. The question is which
ones, and the answer depends on what you want the model to be able to *do*.

A :doc:`../functions/tasks` says exactly that. "Regenerate ATP from glucose
aerobically" or "synthesise L-glutamate from glucose and ammonium" are
statements the model either satisfies or does not, and a task list is the
specification a gap-filled model is held to.

The difficulty is that "the model must perform this task" is difficult to
express as a constraint on an extraction. The trick is to turn it into one
that is not.

The problem
-----------

**Minimise** the number of reactions taken from the universal model,

**subject to**

* **the draft is kept** — every reaction it already has stays, and keeping
  them costs nothing, so only additions are paid for;
* **mass balance** — at steady state nothing piles up;
* **a reaction left out carries no flux** — a reaction not added cannot be
  used by any task;
* **every task is performed** — each one checked in its own conditions, with
  its inputs and outputs opened and its equations added.

The last two matter together: a task is checked with the reactions the model
*has*, so the solver cannot satisfy one task using a reaction it declined to
pay for.

Because only additions are paid for, the answer is "the fewest reactions
that fix the draft", not "the smallest model that performs the tasks" —
those are different problems, and the second would be free to throw away
reactions the draft already had.

Solving this directly means tracking one set of fluxes per task over one
shared set of on/off switches, which is a large problem. The rest of this
page is how SPECTRA avoids paying for it.

The trick: tasks become core reactions
--------------------------------------

For each task, find the reactions without which it becomes infeasible. Their
union is a set of reactions the model must contain, and *that* is something
:func:`~spectra_cobra.spectra_me` already enforces: it is a core set.

That is what :func:`~spectra_cobra.essential_reactions_for_tasks` returns, and
the worked example below runs it on a real model and task list.

The idea of finding a task's essential reactions first and handing them to the
extraction as a core set is taken from ftINIT, where it is what makes
task-constrained extraction affordable enough to run once per cell across a
single-cell dataset:

   Gustafsson, J., Anton, M., Roshanzamir, F., Jörnsten, R., Kerkhoven, E. J.,
   Robinson, J. L., and Nielsen, J. (2023). Generation and analysis of
   context-specific genome-scale metabolic models derived from single-cell
   RNA-Seq data. *Proceedings of the National Academy of Sciences*, 120(6),
   e2217868120. https://doi.org/10.1073/pnas.2217868120

.. warning::

   This is a necessary condition, not a sufficient one. A task with two
   alternative routes has **no** essential reactions at all, because either
   route can be dropped on its own, so forcing the essential set cannot
   guarantee that task survives. It does most of the work for very little
   effort; what tells you whether anything is left over is re-checking the
   tasks afterwards. The last section does exactly that.

Defining tasks
--------------

A task names the metabolites it may take up, the metabolites it may produce,
and any reactions it adds for the duration of the check:

.. code-block:: python

   from spectra_cobra import MetabolicTask, TaskEquation

   ATP_HYDROLYSIS = TaskEquation(
       {"atp_c": -1, "h2o_c": -1, "adp_c": 1, "pi_c": 1, "h_c": 1},
       (1.0, 1000.0),
   )

   task = MetabolicTask(
       "atp_aerobic",
       "Rephosphorylate ATP from glucose, with oxygen",
       inputs={"glc__D_e": (0, 1000), "o2_e": (0, 1000), "pi_e": (0, 1000),
               "h2o_e": (0, 1000), "h_e": (0, 1000)},
       outputs={"co2_e": (0, 1000), "h2o_e": (0, 1000), "h_e": (0, 1000)},
       equations=(ATP_HYDROLYSIS,),
   )

Checking a task closes the model's own boundary, opens exchanges for the
inputs and outputs, adds the equations, and asks whether the result is
feasible. The model itself is never modified.

.. note::

   **A task needs forcing to mean anything.** Bounds default to ``(0, 1000)``
   — allowed, never required — so a task with no positive lower bound
   anywhere is satisfied by doing nothing at all, and passes trivially. Above
   it is the hydrolysis equation's lower bound of 1 that makes the model
   demonstrate it can actually regenerate ATP. For a biosynthesis task, put
   the lower bound on the output instead. :meth:`MetabolicTask.is_forced
   <spectra_cobra.MetabolicTask.is_forced>` reports whether anything forces
   a task.

Mark the negative controls ``should_fail=True``. A model that can make ATP
from nothing is leaking, and the control is how you find out.

A worked example
----------------

On the textbook model, with four tasks:

.. code-block:: python

   from cobra.io import load_model
   from spectra_cobra import (
       MetabolicTask, TaskEquation, check_tasks,
       essential_reactions_for_tasks, spectra_me,
   )

   universal = load_model("textbook")

   HYDROLYSIS = {"atp_c": -1.0, "h2o_c": -1.0,
                 "adp_c": 1.0, "pi_c": 1.0, "h_c": 1.0}
   AEROBIC = ["glc__D_e", "o2_e", "nh4_e", "pi_e", "h2o_e", "h_e"]
   RESPIRED = {m: (0.0, 1000.0) for m in
               ["co2_e", "h2o_e", "h_e", "ac_e", "etoh_e",
                "for_e", "lac__D_e", "succ_e"]}
   FERMENTED = {m: (0.0, 1000.0) for m in
                ["co2_e", "h2o_e", "h_e", "ac_e", "etoh_e", "for_e"]}

   tasks = [
       MetabolicTask(
           id="atp_aerobic",
           inputs={m: (0.0, 10.0) for m in AEROBIC},
           outputs=RESPIRED,
           equations=[TaskEquation(HYDROLYSIS, (1.0, 1000.0))],
       ),
       MetabolicTask(
           id="atp_anaerobic",
           inputs={m: (0.0, 10.0) for m in AEROBIC if m != "o2_e"},
           outputs=FERMENTED,
           equations=[TaskEquation(HYDROLYSIS, (1.0, 1000.0))],
       ),
       MetabolicTask(
           id="make_glutamate",
           inputs={m: (0.0, 10.0) for m in AEROBIC},
           outputs=RESPIRED,
           equations=[TaskEquation({"glu__L_c": -1.0}, (1.0, 1000.0))],
       ),
       # A negative control: nothing may be taken up, so ATP cannot be made.
       MetabolicTask(
           id="no_free_atp",
           inputs={},
           outputs=RESPIRED,
           equations=[TaskEquation(HYDROLYSIS, (1.0, 1000.0))],
           should_fail=True,
       ),
   ]

   for result in check_tasks(universal, tasks):
       print(result.task.id, result.feasible, result.ok)

.. code-block:: text

   atp_aerobic      feasible=True  as expected=True
   atp_anaerobic    feasible=True  as expected=True
   make_glutamate   feasible=True  as expected=True
   no_free_atp      feasible=False as expected=True

Always do this before gap-filling. A task the universal model cannot perform
can never be gap-filled towards, because there is nothing to add that would
help, and it will sit in your results as a permanent failure.

Now the essential reactions:

.. code-block:: python

   core, per_task = essential_reactions_for_tasks(universal, tasks)

   for task_id, reactions in sorted(per_task.items()):
       print(f"{task_id:<16} {len(reactions):>2} essential")
   print(f"union (the core set): {len(core)}")

.. code-block:: text

   atp_aerobic       1 essential
   atp_anaerobic    12 essential
   make_glutamate   10 essential
   union (the core set): 17

Aerobic ATP production has exactly one essential reaction, because almost
every step of it has an alternative; anaerobic ATP has twelve, because
fermentation is a narrow path. That spread is the warning above, in numbers.

Knock four reactions out to make a draft:

.. code-block:: python

   draft = universal.copy()
   draft.remove_reactions(["CYTBD", "ACALD", "GLUSy", "AKGDH"],
                          remove_orphans=False)

   failing = [result.task.id for result in check_tasks(draft, tasks)
              if not result.ok]
   print(f"tasks the draft fails: {failing}")

.. code-block:: text

   tasks the draft fails: ['atp_anaerobic', 'make_glutamate']

Two of the three break, not all three. ``atp_aerobic`` survives, because its
one essential reaction was not among the four removed and the rest of aerobic
ATP production has alternatives to spare. That is the warning above happening
in practice: the essential set is a necessary condition, and a task with
alternative routes contributes few essential reactions, or none.

Gap-filling
-----------

:func:`~spectra_cobra.gapfill_for_tasks` does the whole thing: find the
essential reactions, extract around them, and apply the additions on top of
the draft.

.. code-block:: python

   from spectra_cobra import gapfill_for_tasks

   result = gapfill_for_tasks(draft, universal, tasks, tol=1e-3, seed=0)

   print(f"added:       {sorted(result.added)}")
   print(f"satisfied:   {sorted(result.satisfied)}")
   print(f"unsatisfied: {sorted(result.unsatisfied)}")

.. code-block:: text

   added:       ['ACALD']
   satisfied:   ['atp_aerobic', 'atp_anaerobic', 'make_glutamate']
   unsatisfied: []

One reaction was enough. ``CYTBD``, ``GLUSy`` and ``AKGDH`` were not put back,
because the network has other ways to respire, to make glutamate and to run
the TCA cycle. That is the point of gap-filling towards tasks rather than
towards a reaction list: it restores the *ability*, not the reactions you
happened to delete.

The result carries the model as well as the account, so you can check it
independently rather than taking ``satisfied`` on trust:

.. code-block:: python

   still_failing = [r.task.id for r in check_tasks(result.model, tasks)
                    if not r.ok]
   print(f"tasks failing now: {still_failing}")
   print(f"{len(result.model.reactions)} reactions, draft had {len(draft.reactions)}")

.. code-block:: text

   tasks failing now: []
   92 reactions, draft had 91

Arguments worth knowing
-----------------------

``problem_type``
   The formulation the **core solve** uses -- the one
   :func:`~spectra_cobra.spectra_me` runs over the universal model. It is
   ``"minNetMILP"`` by default here, unlike
   :func:`~spectra_cobra.spectra_me` itself: gap-filling is usually run once
   on one draft, and the exact answer is worth the mixed-integer solve. Pass
   ``"minNetLP"`` if it is too slow.

   It does **not** govern the repair below, which is always mixed-integer
   whatever you pass here. Repairing one failing task is a small problem, so
   there is nothing to gain from relaxing it.

``repair``
   Whether to re-check the tasks afterwards and add what is still missing
   (default True). The core set is a necessary condition, not a sufficient
   one, so an extraction around it can leave a task unsatisfied; repair is
   what closes that gap. It costs one task check plus, if anything is still
   failing, a second solve. On this example it changes nothing -- the core
   set was sufficient -- but it is what makes the ``satisfied`` list
   trustworthy in general.

   The repair solve is **always** ``minNetMILP``, whatever ``problem_type``
   you passed for the core solve. Repairing one failing task is a small
   problem -- it runs inside that task's own constraints, over what the
   draft is still missing -- so the exact answer is cheap, and there is
   nothing to gain by relaxing it to an LP.

``essential_reactions``
   Pass a core set you have already computed, to avoid paying for it twice.

``free_reactions``
   Reactions to charge nothing for. The draft's own reactions are free
   already; this is for anything else you want the extraction not to count.

.. important::

   :func:`~spectra_cobra.spectra_me` on its own is an *extraction*: it returns
   a minimal network around the core and discards everything else, including
   most of your draft. Gap-filling is the extraction's additions applied **on
   top of** the draft, and :func:`~spectra_cobra.gapfill_for_tasks` is what
   does that union for you. Doing it by hand is easy to get subtly wrong.

At genome scale
---------------

The same workflow on iJO1366 (2583 reactions), with 23 tasks — aerobic and
anaerobic ATP production, growth on acetate, and one biosynthesis task per
proteinogenic amino acid:

.. code-block:: python

   AMINO_ACIDS = ["ala__L", "arg__L", "asn__L", "asp__L", "cys__L", "gln__L",
                  "glu__L", "gly", "his__L", "ile__L", "leu__L", "lys__L",
                  "met__L", "phe__L", "pro__L", "ser__L", "thr__L", "trp__L",
                  "tyr__L", "val__L"]
   AEROBIC = ["glc__D_e", "o2_e", "nh4_e", "pi_e", "h2o_e", "h_e",
              "so4_e", "k_e"]

   universal = load_model("iJO1366")
   tasks = [
       MetabolicTask(id="atp_aerobic",
                     inputs={m: (0.0, 10.0) for m in AEROBIC}, outputs=RESPIRED,
                     equations=[TaskEquation(HYDROLYSIS, (1.0, 1000.0))]),
       MetabolicTask(id="atp_anaerobic",
                     inputs={m: (0.0, 10.0) for m in AEROBIC if m != "o2_e"},
                     outputs=RESPIRED,
                     equations=[TaskEquation(HYDROLYSIS, (1.0, 1000.0))]),
       MetabolicTask(id="grow_on_acetate",
                     inputs={**{m: (0.0, 10.0) for m in AEROBIC
                                if m != "glc__D_e"}, "ac_e": (0.0, 10.0)},
                     outputs=RESPIRED,
                     equations=[TaskEquation(HYDROLYSIS, (1.0, 1000.0))]),
   ] + [
       MetabolicTask(id=f"make_{aa}",
                     inputs={m: (0.0, 10.0) for m in AEROBIC}, outputs=RESPIRED,
                     equations=[TaskEquation({f"{aa}_c": -1.0}, (0.1, 1000.0))])
       for aa in AMINO_ACIDS
   ]

All 23 are feasible on iJO1366 itself, and they have 95 essential reactions
between them. The draft deletes 150 random non-exchange reactions, two of
which were essential to some task:

.. code-block:: python

   import random

   internal = sorted(r.id for r in universal.reactions
                     if r not in universal.boundary)
   draft = universal.copy()
   draft.remove_reactions(random.Random(0).sample(internal, 150),
                          remove_orphans=False)

Checking the tasks took 0.4 s and finding the 95 essential reactions 5.1 s,
so the whole preparation is cheaper than a single genome-scale MILP.

Two of the 23 tasks are skipped rather than attempted: ``make_cys__L`` and
``make_met__L`` are infeasible on iJO1366 itself, so no gap-fill could help
them. The last section of this page explains why.

Which objective?
----------------

Those three rows are one function call each, differing only in arguments:

.. code-block:: python

   from spectra_cobra import gapfill_for_tasks

   runs = {
       "tradeOff": dict(problem_type="tradeOff"),
       "minNetMILP, repair=False": dict(problem_type="minNetMILP", repair=False),
       "minNetMILP + repair": dict(problem_type="minNetMILP", repair=True),
   }

   for name, options in runs.items():
       result = gapfill_for_tasks(draft, universal, tasks, **options)
       print(f"{name:<26} {len(result.model.reactions):>5} reactions, "
             f"{len(result.added):>4} added, "
             f"{len(result.unsatisfied)} failing")

.. code-block:: text

   tradeOff                    2522 reactions,   89 added, 0 failing
   minNetMILP, repair=False    2448 reactions,   15 added, 1 failing
   minNetMILP + repair         2449 reactions,   16 added, 0 failing

.. list-table::
   :header-rows: 1
   :widths: 34 16 16 16 18

   * - Model
     - Reactions
     - Added
     - Tasks failing
     - Time
   * - draft
     - 2433
     - —
     - 6
     - —
   * - ``tradeOff``
     - 2522
     - 89
     - **0**
     - 7.7 s
   * - ``minNetMILP``, ``repair=False``
     - 2448
     - 15
     - 1
     - 6.4 s
   * - ``minNetMILP`` + ``repair``
     - 2449
     - 16
     - **0**
     - 7.8 s

The draft row is the starting point, and the times are wall clock on one
machine — useful for their ratio rather than their absolute value.

The first two differ by a factor of seven in what they add, and the
difference is instructive rather than a matter of taste.

``minNetMILP`` adds the minimum: 16 reactions, essentially just the deleted
essential ones. One task — using acetate as a carbon source — still fails,
because it has alternative routes and therefore contributed nothing to the
core set. This is the warning at the top of the page, happening.

``tradeOff`` adds 109 and fixes everything, but not because it understood the
problem better. Signed weights make it *reward* including draft reactions, so
it builds a larger active network, and the extra reactions happen to restore
the acetate route. That is luck, and on another draft it would not be.

``tradeOff`` works here, but note the condition that makes it work: the core
is a set of metabolic reactions. Make biomass the core instead and
``tradeOff`` becomes infeasible outright, because it requires every included
reaction to carry at least ``tol`` while growth needs trace fluxes far below
that — on iJO1366 at full growth the smallest necessary flux is around
``2e-06``. Use ``minNetMILP`` whenever growth is what you are asking for; its
binaries pin a reaction to zero when it is out but impose no floor when it is
in, so a trace flux is perfectly acceptable.

``repair`` is what closes the gap the core set could not express. It re-checks
the tasks after the extraction and, for any still failing, solves once more
inside that task's own constraints — where the model is already forced to
perform it, so minimising the reactions not already present gives the
cheapest repair directly. Here it cost one extra reaction and one extra
solve, and took the failures to zero: 2449 reactions against ``tradeOff``'s
2522. Minimise, then repair what the core could not express — 73 reactions
smaller, and correct for a reason rather than by accident.

Where this fits
---------------

The two phases serve two different jobs, and only one of them should call
:func:`~spectra_cobra.gapfill_for_tasks` end to end.

**Automated genome-scale reconstruction.** A draft comes out of a pipeline
and a task list says what the organism must be able to do. Call
:func:`~spectra_cobra.gapfill_for_tasks` directly, as above: it derives the
essential reactions, extracts around them, and repairs what is left.

**Context-specific modelling.** Here the core set is not the tasks alone. It
is the reactions your evidence calls active *plus* the task-essential ones —
what the data says is on, together with what the model must keep to stay
functional. Do not hand this to
:func:`~spectra_cobra.gapfill_for_tasks`: its core solve would recompute the
essential reactions and extract again with flat weights, discarding the
evidence that made the model context-specific in the first place. Run the
extraction yourself, then repair:

.. code-block:: python

   from spectra_cobra import (
       essential_reactions_for_tasks, gapfill_for_tasks, spectra_me,
   )

   # 1. What the tasks require, whatever the data says.
   essential, _ = essential_reactions_for_tasks(universal, tasks)

   # 2. What your context says is active, and how strongly.
   expressed = {"PGI", "PFK", "FBA", "TPI", "GAPD", "PGK"}
   evidence = {rxn.id: (0.1 if rxn.id in expressed else 1.0)
               for rxn in universal.reactions}

   # 3. One core set, one extraction, your evidence as the weights.
   core = sorted(expressed | essential)
   context = spectra_me(universal, core, tol=1e-4, weights=evidence, seed=0)
   print(f"context-specific model: {len(context.reactions)} reactions")

   # 4. Repair only -- the core solve has already happened.
   result = gapfill_for_tasks(context, universal, tasks,
                              essential_reactions=essential,
                              core_solve=False)
   print(f"after repair: {len(result.model.reactions)} reactions, "
         f"{len(result.satisfied)}/{len(tasks)} tasks")

.. code-block:: text

   context-specific model: 207 reactions
   after repair: 219 reactions, 23/23 tasks

``core_solve=False`` is what stops the second call redoing the first. The
repair still runs, so the tasks the core set could not guarantee are still
caught — which matters more here than in the automated case, because a core
set built from expression evidence was never chosen to satisfy the tasks.

Validate the task list first
----------------------------

Two of the first 23 tasks tried here were infeasible on iJO1366 itself. The
cause was specific: both asked for the *extracellular* amino acid, and
iJO1366 synthesises glutamine and methionine but has no route to export
them. ``GLNabcpp`` and ``METabcpp`` are irreversible inward, so the maximum
production of ``gln__L_e`` is exactly zero while ``gln__L_c`` reaches 11.5.

No gap-fill can repair that, because nothing in the universal model would
help. Writing the tasks against the cytosolic metabolite, which is what
"can this network synthesise X" actually means, made all 23 feasible.

Next
----

* :doc:`gapfilling_media` — the same idea, with growth in each of a panel of
  media as the requirement instead of a task list
* :doc:`formulations` — what the objectives optimise
* :doc:`../functions/tasks` — the task API
