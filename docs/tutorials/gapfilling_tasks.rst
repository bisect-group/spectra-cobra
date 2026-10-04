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

The difficulty is that "the model must perform this task" is awkward to
express as a constraint on an extraction. The trick is to turn it into one
that is not.

The trick: tasks become core reactions
--------------------------------------

For each task, find the reactions without which it becomes infeasible. Their
union is a set of reactions the model must contain, and *that* is something
:func:`~spectra_cobra.spectra_me` already enforces: it is a core set.

.. code-block:: python

   from spectra_cobra import essential_reactions_for_tasks

   core, per_task = essential_reactions_for_tasks(universal, tasks)

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
   from spectra_cobra import check_tasks, essential_reactions_for_tasks, spectra_me

   universal = load_model("textbook")

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

.. code-block:: text

   atp_aerobic       1 essential
   atp_anaerobic    12 essential
   make_glutamate   10 essential
   union (the core set): 17

Aerobic ATP production has exactly one essential reaction, because almost
every step of it has an alternative; anaerobic ATP has twelve, because
fermentation is a narrow path. That spread is the warning above, in numbers.

Knock four reactions out to make a draft, and all three tasks break:

.. code-block:: python

   draft = universal.copy()
   draft.remove_reactions(["CYTBD", "ACALD", "GLUSy", "AKGDH"],
                          remove_orphans=False)

.. code-block:: text

   tasks the draft fails: ['atp_aerobic', 'atp_anaerobic', 'make_glutamate']

Gap-fill by extracting around the core, rewarding the reactions the draft
already has and charging for everything else:

.. code-block:: python

   draft_ids = {rxn.id for rxn in draft.reactions}
   weights = {rxn.id: (1.0 if rxn.id in draft_ids else -1.0)
              for rxn in universal.reactions}

   filled = spectra_me(universal, sorted(core), tol=1e-3, weights=weights,
                       problem_type="tradeOff", seed=0)

   # The extraction returns a minimal network around the core. The gap-filled
   # model is that network unioned with the draft: the point is to add to the
   # draft, not to replace it.
   added = {rxn.id for rxn in filled.reactions} - draft_ids

.. code-block:: text

   added back: ['ACALD', 'CYTBD']
   tasks failing now: []

Two of the four removed reactions were enough. ``GLUSy`` and ``AKGDH`` were
not put back because the network has other ways to make glutamate and to run
the TCA cycle, which is the point of gap-filling towards tasks rather than
towards a particular reaction list.

.. important::

   The union with the draft is not optional. :func:`~spectra_cobra.spectra_me`
   is an *extraction*: on its own it returns a minimal network around the core
   and discards everything else, including most of your draft. Gap-filling is
   the extraction's additions applied **on top of** the draft.

At genome scale
---------------

The same workflow on iJO1366 (2583 reactions), with 23 tasks — energy
metabolism, acetate utilisation, and one biosynthesis task per proteinogenic
amino acid — and a draft made by deleting 150 random non-exchange reactions,
4 of which were essential to some task:

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
     - 5
     - —
   * - ``tradeOff``
     - 2476
     - 43
     - **0**
     - 6.1 s
   * - ``minNetMILP``
     - 2441
     - 8
     - 1
     - 4.4 s
   * - ``minNetMILP`` + repair
     - 2442
     - 9
     - **0**
     - +1 solve

Checking the tasks took 0.5 s and finding the 93 essential reactions 5.3 s,
so the whole preparation is cheaper than a single genome-scale MILP.

Which objective?
----------------

The two rows above differ by a factor of five in what they add, and the
difference is instructive rather than a matter of taste.

``minNetMILP`` adds the minimum: 8 reactions, essentially just the deleted
essential ones. One task — using acetate as a carbon source — still fails,
because it has alternative routes and therefore contributed nothing to the
core set. This is the warning at the top of the page, happening.

``tradeOff`` adds 43 and fixes everything, but not because it understood the
problem better. Signed weights make it *reward* including draft reactions, so
it builds a larger active network, and the extra reactions happen to restore
the acetate route. That is luck, and on another draft it would not be.

``tradeOff`` works here, but note the condition that makes it work: the core
is a set of metabolic reactions. Make biomass the core instead and
``tradeOff`` becomes infeasible outright, because it requires every included
reaction to carry at least ``tol`` and growth needs trace fluxes far below
that. :doc:`gapfilling_media` has the measurement.

So prefer ``minNetMILP`` and then repair what is left:

.. code-block:: python

   from spectra_cobra.formulations import min_net_milp
   from spectra_cobra.tasks import task_constraints

   for result in check_tasks(gapfilled, tasks):
       if result.ok:
           continue
       # Inside the task's own constraints the model is already forced to
       # perform it, so minimising the count of reactions that are not
       # already present gives the cheapest repair directly.
       with task_constraints(universal, result.task) as constrained:
           weights = {rxn.id: (0.0 if rxn.id in present else 1.0)
                      for rxn in constrained.reactions}
           kept = min_net_milp(constrained, {}, weights, 1e-4, True, 300.0, None)

On the run above that added **one** reaction and brought the failures to
zero, for a final model of 2442 reactions against ``tradeOff``'s 2476. Minimise,
then repair what the core could not express: 34 reactions smaller, and correct
for a reason rather than by accident.

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

* :doc:`gapfilling_media` — the same idea, with growth in a medium as the
  requirement instead of a task list
* :doc:`formulations` — what the objectives optimise
* :doc:`../functions/tasks` — the task API
