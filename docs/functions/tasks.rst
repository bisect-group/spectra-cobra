Metabolic tasks
===============

A *task* is a statement about what a network must be able to do: take up
these metabolites, produce those, carry flux through this reaction. It is
checked by closing the model's boundary, opening exchanges for the task's
inputs and outputs, adding the task's own reactions, and asking whether the
result is feasible.

Usage
-----

.. code-block:: python

   from spectra_cobra import MetabolicTask, TaskEquation, check_tasks

   task = MetabolicTask(
       "atp_aerobic",
       "Rephosphorylate ATP from glucose, with oxygen",
       inputs={"glc__D_e": (0, 1000), "o2_e": (0, 1000), "pi_e": (0, 1000),
               "h2o_e": (0, 1000), "h_e": (0, 1000)},
       outputs={"co2_e": (0, 1000), "h2o_e": (0, 1000), "h_e": (0, 1000)},
       equations=(TaskEquation(
           {"atp_c": -1, "h2o_c": -1, "adp_c": 1, "pi_c": 1, "h_c": 1},
           (1.0, 1000.0)),),
   )

   for result in check_tasks(model, [task]):
       print(result.task.id, result.feasible, result.ok)

Bounds default to ``(0, 1000)`` — allowed, never required — so something must
carry a positive lower bound or the task is satisfied by the all-zero flux
vector. Above that is the equation; for a biosynthesis task it is usually the
output.

Why this is here
----------------

Requiring a model to *perform a task* is awkward to express as a constraint
on an extraction. Requiring it to *contain a reaction* is exactly what
:doc:`spectra_me` already does. So compute, for each task, the reactions
without which it becomes infeasible, and hand their union to the extraction
as the core set:

.. code-block:: python

   from spectra_cobra import essential_reactions_for_tasks

   core, per_task = essential_reactions_for_tasks(universal, tasks)
   extracted = spectra_me(universal, sorted(core), tol=1e-4)

This is a necessary condition rather than a sufficient one: a task with two
alternative routes has no essential reactions, so forcing the set cannot
guarantee the task survives. See :doc:`../tutorials/gapfilling_tasks` for
what to do about the remainder.

Cost
----

Only a reaction carrying flux in a feasible solution can be essential, since
a solution that leaves it at zero stays available when it is fixed to zero.
The search therefore starts from the reactions one solution uses and shrinks
that set by repeatedly minimising the total flux through whatever is left,
which rules out in a single LP every reaction some other route can avoid.
Only the survivors are tested individually. On iJO1366, 23 tasks took 5.3 s.

Reference
---------

.. autoclass:: spectra_cobra.MetabolicTask
   :members:

.. autoclass:: spectra_cobra.TaskEquation

.. autoclass:: spectra_cobra.TaskResult

.. autofunction:: spectra_cobra.check_task

.. autofunction:: spectra_cobra.check_tasks

.. autofunction:: spectra_cobra.essential_reactions_for_task

.. autofunction:: spectra_cobra.essential_reactions_for_tasks

See also
--------

* :doc:`../tutorials/gapfilling_tasks` — gap-filling towards a task list
* :doc:`../tutorials/gapfilling_media` — gap-filling towards growth in a medium
* :doc:`spectra_me` — the extraction the core set feeds
