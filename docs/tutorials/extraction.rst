Extracting a context-specific model
===================================

Given a universal model and a set of reactions you believe are active in your
context, :func:`~spectra_cobra.spectra_me` returns the smallest surrounding
network that lets them all carry flux.

The basic case
--------------

.. code-block:: python

   from cobra.io import load_model
   from spectra_cobra import spectra_cc, spectra_me

   model = load_model("textbook")
   consistent = spectra_cc(model, tol=1e-3)

   extracted = spectra_me(consistent, core_reactions=["PGI"], tol=1e-3, seed=0)
   print(len(extracted.reactions))
   print("PGI" in {r.id for r in extracted.reactions})

.. code-block:: text

   27
   True

From 87 consistent reactions down to 27, and the core reaction is in there.
Note the order: consistency first, extraction second.
:func:`~spectra_cobra.spectra_me` requires a consistent model, because it has
to be able to drive every core reaction past ``tol`` simultaneously.

The core set grows the model as needed
--------------------------------------

.. code-block:: python

   for core in (["PGI"], ["PGI", "PFK", "FBA"], ["SUCDi", "FUM"]):
       extracted = spectra_me(consistent, core, tol=1e-3, seed=0)
       print(core, "->", len(extracted.reactions), "reactions")

.. code-block:: text

   ['PGI'] -> 27 reactions
   ['PGI', 'PFK', 'FBA'] -> 27 reactions
   ['SUCDi', 'FUM'] -> 37 reactions

Adding ``PFK`` and ``FBA`` to ``PGI`` costs nothing — they are on the same
glycolytic route, which was already pulled in. Asking for two TCA reactions
instead costs ten more, since that part of the network needs more support.

Specifying the core set
-----------------------

Identifiers, :class:`cobra.Reaction` objects and integer indices all work, so
you can hand over whatever your upstream analysis produced:

.. code-block:: python

   spectra_me(consistent, ["PGI"], tol=1e-3)                       # identifier
   spectra_me(consistent, [consistent.reactions.PGI], tol=1e-3)    # object
   spectra_me(consistent, [12], tol=1e-3)                          # index

A reaction that is not in the model raises
:class:`~spectra_cobra.SpectraError` rather than being skipped quietly.

When the core set is impossible
-------------------------------

If no flux distribution lets every core reaction carry ``tol`` at once, you
get an error naming the problem:

.. code-block:: python

   from spectra_cobra import SpectraInfeasibleCoreError

   # 'FRUpts2' is blocked in the textbook model's default medium.
   try:
       spectra_me(model, ["PGI", "FRUpts2"], tol=1e-3)
   except SpectraInfeasibleCoreError as error:
       print(error)

.. code-block:: text

   No flux distribution lets all 1 pinned irreversible core reactions carry
   at least tol=0.001 at once (solver status 'infeasible'). Either the core
   set is mutually inconsistent, or some of its reactions are blocked in this
   model; spectra_ccme handles the latter by dropping them.

Two causes, as the message says. Either a core reaction is blocked in the
model — which is what happens if you skipped the consistency check — or the
core reactions are individually fine but mutually incompatible. The count is
of *pinned* reactions, which is why it says one rather than two: ``PGI`` is
reversible and so handled through the objective, while only the irreversible
``FRUpts2`` is pinned through its bounds.

Handling an inconsistent universal model
----------------------------------------

Rather than checking consistency yourself and intersecting the result with
your core set, let :func:`~spectra_cobra.spectra_ccme` do both:

.. code-block:: python

   from spectra_cobra import spectra_ccme

   extracted, blocked_core = spectra_ccme(model, ["PGI", "FRUpts2"], tol=1e-3,
                                          seed=0)
   print(f"{len(extracted.reactions)} reactions")
   print(f"dropped as blocked: {blocked_core}")

.. code-block:: text

   27 reactions
   dropped as blocked: ['FRUpts2']

It reports what it dropped rather than failing, which is usually what you
want when the core set came from omics data and may name reactions the
universal model cannot support.

Weights
-------

Weights steer which reactions get pulled in. They are a dict keyed by reaction
identifier, and anything you leave out defaults to 1.0:

.. code-block:: python

   weights = {"SUCDi": 100.0, "FUM": 100.0}   # expensive, avoid if possible
   extracted = spectra_me(consistent, ["PGI"], tol=1e-3, weights=weights)

What a weight *means* depends on which formulation you use, and the sign
convention is not the same for all of them — see :doc:`formulations`.

Dropping orphaned genes
-----------------------

Removing reactions can leave genes with nothing to act on. They are dropped
by default; ``remove_genes=False`` keeps them:

.. code-block:: python

   extracted = spectra_me(consistent, ["PGI"], tol=1e-3, remove_genes=False)

Keeping them is the slower of the two: the extracted model is built from the
reactions kept, so the orphans have to be put back afterwards.

Reproducibility
---------------

The direction phase breaks ties between core reactions with randomised
objective coefficients, so an unseeded run returns a different model each
time. ``seed`` pins them, and it does so across processes as well as within
one:

.. code-block:: python

   core = ["PGI", "PFK", "FBA"]

   a = spectra_me(consistent, core, tol=1e-3, seed=0)
   b = spectra_me(consistent, core, tol=1e-3, seed=0)
   assert {r.id for r in a.reactions} == {r.id for r in b.reactions}

Core identifiers are sorted before the LP is built, so each reaction keeps
its own coefficient whatever ``PYTHONHASHSEED`` does to the iteration order
of the core set. Without that, a seeded run reproduces only inside a single
process. The solver is the remaining source of variation.

Which formulation to extract with
---------------------------------

Everything above used the default, ``problem_type="minNetLP"``: one LP,
minimising total absolute flux. ``"minNetMILP"`` asks the same question
exactly rather than through a relaxation and so returns a smaller model --
16 reactions against 27 for the core set used here -- but it is a
mixed-integer solve, and that is where the cost lands.

All four formulations keep every core reaction. The core is pinned by the
directions settled in phase one, not by the objective, so the choice changes
what is kept *around* the core and never whether the core survives. What it
also changes is the meaning of a weight, which is the part most worth
reading up on before switching: see :doc:`formulations`.

The size of the model decides whether the exact answer is affordable. On the
textbook model every formulation finishes in well under a tenth of a second
and you may as well use ``minNetMILP``. On a genome-scale model it does not
finish at all without a ``time_limit`` -- see :doc:`recon3d`, where it is
given ten minutes against an LP that takes seconds.

A useful rule, then, is to let the number of models decide:

* **Many models from one universal model** -- a context-specific model per
  cell in a single-cell dataset, per sample, or per condition -- favours the
  **LP formulations**. A per-model cost that is tolerable once becomes the
  whole runtime when it is paid thousands of times, and a slightly larger
  model is rarely what limits that kind of analysis.
* **One model you intend to study closely** favours **minNetMILP**, with a
  ``time_limit`` set if the model is large. Here the size of the result is
  the result, and the minutes are paid once.

Next
----

* :doc:`formulations` — what to optimise, and what weights mean
* :doc:`alternative_solutions` — several models instead of one
* :doc:`recon3d` — the same workflow at genome scale
