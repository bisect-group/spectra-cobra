Alternative solutions
=====================

A core set rarely has one answer. Several different networks can support the
same reactions equally well, and which one you get depends on arbitrary
tie-breaking inside the solver. ``n_solutions`` asks for more than one.

.. code-block:: python

   models = spectra_me(consistent, ["PGI"], tol=1e-3, n_solutions=3)

With ``n_solutions=1`` (the default) you get a model. With more, you get a
list of models.

Two methods
-----------

``"pathwayExclusion"``
   Adds a constraint forbidding each reaction set already returned, so every
   model is genuinely different. Needs a mixed-integer ``problem_type``.

``"coreDirection"`` (the default)
   Re-runs the direction phase with the two LPs in a random order, so the core
   reactions may run the other way and pull in a different network. Works with
   any ``problem_type``, but the alternatives are not guaranteed to differ.

Pathway exclusion
-----------------

.. code-block:: python

   models = spectra_me(consistent, ["PGI"], tol=1e-3, n_solutions=3,
                       alt_solution_method="pathwayExclusion",
                       problem_type="minNetMILP", seed=0)

   print(len(models), [len(m.reactions) for m in models])
   print(len({frozenset(r.id for r in m.reactions) for m in models}))

.. code-block:: text

   3 [16, 18, 19]
   3

Three models, three distinct reaction sets, and they grow in size: the first
is the true minimum at 16, and each later one is the best available once the
earlier answers are ruled out.

.. important::

   **You may get fewer models than you asked for.** Once the exclusions leave
   no feasible answer, the network has no further alternative to offer, and
   you get the ones that exist rather than an error. Always check ``len()``
   on the result rather than assuming:

   .. code-block:: python

      models = spectra_me(consistent, ["PGI"], tol=1e-3, n_solutions=50,
                          alt_solution_method="pathwayExclusion",
                          problem_type="minNetMILP", seed=0)
      print(f"asked for 50, got {len(models)}")

   .. code-block:: text

      asked for 50, got 50

   The textbook model has far more than fifty ways to support ``PGI``, so
   that one does not run out. A small network does: on the three-pathway toy
   model, asking for five returns exactly the three routes that exist.

Excluding solutions you already have
------------------------------------

``previous_solutions`` seeds the exclusion set, which is useful for resuming
an enumeration or for ruling out a model from an earlier run:

.. code-block:: python

   first = spectra_me(consistent, ["PGI"], tol=1e-3,
                      problem_type="minNetMILP", seed=0)
   already = {r.id for r in first.reactions}

   nxt = spectra_me(consistent, ["PGI"], tol=1e-3,
                    problem_type="minNetMILP",
                    previous_solutions=[already], seed=0)

   assert {r.id for r in nxt.reactions} != already

Passing it with an LP formulation raises
:class:`~spectra_cobra.SpectraError`, because only the mixed-integer
formulations have the binary variables the constraint is written over.

Core direction
--------------

.. code-block:: python

   models = spectra_me(consistent, ["PGI"], tol=1e-3, n_solutions=3, seed=0)
   print([len(m.reactions) for m in models])

   distinct = {frozenset(r.id for r in m.reactions) for m in models}
   print(f"{len(models)} models, {len(distinct)} distinct")

.. code-block:: text

   [27, 36, 27]
   3 models, 2 distinct

This is the cheaper option — it needs no MILP — but as that output shows, it
varies the core directions rather than forbidding solutions outright, so two
of the three landed on the same network. De-duplicate if it matters to you,
or use pathway exclusion.

Reproducibility
---------------

Both methods are randomised, so pass ``seed`` for a reproducible set:

.. code-block:: python

   a = spectra_me(consistent, ["PGI"], tol=1e-3, n_solutions=3, seed=42)
   b = spectra_me(consistent, ["PGI"], tol=1e-3, n_solutions=3, seed=42)

   assert [{r.id for r in m.reactions} for m in a] == \
          [{r.id for r in m.reactions} for m in b]

Next
----

* :doc:`formulations` — which ``problem_type`` supports what
* :doc:`recon3d` — a genome-scale walkthrough
