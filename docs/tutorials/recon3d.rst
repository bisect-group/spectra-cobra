Genome scale: Recon3D
=====================

Everything so far used toy models and the 95-reaction textbook model. This
page runs the same workflow on Recon3D, the human reconstruction SPECTRA's own
experiments use, and reports what actually happens — including where the
default formulation falls short.

The model used here is the consistent Recon3D from the MATLAB repository:
11303 reactions over 6388 metabolites.

Consistency at genome scale
---------------------------

.. code-block:: python

   from spectra_cobra import consistent_reaction_ids

   model.solver = "gurobi"
   ids, n_lps = consistent_reaction_ids(model, tol=1e-4, seed=0)
   print(f"{len(ids)}/{len(model.reactions)} consistent in {n_lps} LPs")

.. code-block:: text

   11303/11303 consistent in 8 LPs

All of it, in eight LPs and about fifteen seconds. Two things to take from
that. The model is genuinely flux consistent, which confirms the ``cons``
in its name. And the LP count barely moves with model size — the textbook
model also took eight.

Extraction with a random core set
---------------------------------

A realistic test: draw a random core set, extract, and check the two things
the method promises — that every core reaction is in the result, and that the
result is flux consistent.

.. code-block:: python

   import numpy as np
   from cobra.flux_analysis import find_blocked_reactions
   from spectra_cobra import spectra_me

   rng = np.random.default_rng(0)
   pool = sorted(r.id for r in model.reactions)
   core = sorted(rng.choice(pool, size=200, replace=False).tolist())

   extracted = spectra_me(model, core, tol=1e-4, problem_type="minNetLP",
                          seed=0)

   kept = {r.id for r in extracted.reactions}
   print(f"{len(kept)} reactions")
   print(f"core kept: {sum(c in kept for c in core)}/{len(core)}")
   print(f"blocked: {len(find_blocked_reactions(extracted))}")

Across twelve trials at four core sizes, the results were:

.. list-table::
   :header-rows: 1
   :widths: 12 12 18 18 20 20

   * - core
     - trials
     - extracted
     - core kept
     - blocked
     - of which core
   * - 10
     - 3
     - 167–177
     - all
     - 0, 0, 0
     - —
   * - 50
     - 3
     - 525–927
     - all
     - 0, 0, **112**
     - 4
   * - 200
     - 3
     - 1495–1791
     - all (one failed)
     - **37**, **19**
     - 2, 2
   * - 500
     - 3
     - 2752–2792
     - all
     - **15**, **16**, 0
     - 5, 3

Each run takes roughly a minute, most of it building the LP rather than
solving it.

.. warning::

   **``minNetLP`` does not guarantee a flux consistent model.** Every core
   reaction is always *present* in the result, but in most trials a handful of
   reactions in the extracted model — **including some core reactions** —
   cannot carry flux in it. This package's consistency check and cobrapy's
   independent FVA-based ``find_blocked_reactions`` agree exactly on every
   count, so it is not a detection artefact.

   The cause is the inclusion rule, which comes from MATLAB:
   a reaction is kept if its flux in the extraction LP exceeds
   ``tol * 1e-7``, floored here at ``model.tolerance``, so around
   ``1e-7``. An L1 objective has many optimal solutions and spreads tiny
   fluxes across thousands of reactions, so dropping the ten thousand
   reactions below that threshold discards flux that was doing real balancing
   work. On a hub metabolite such as ``h[c]``, which 2079 reactions touch,
   that discarded dust sums to more than ``tol`` — so the kept flux vector is
   only approximately mass balanced, and the reactions whose flux was itself
   near ``tol`` can no longer carry it.

   **Always check the result** rather than assuming it:

   .. code-block:: python

      from spectra_cobra import blocked_reaction_ids

      blocked = set(blocked_reaction_ids(extracted, tol=1e-4))
      assert not blocked, f"{len(blocked)} blocked, {len(blocked & set(core))} core"

Getting a consistent model
--------------------------

Two options, with different costs.

**Use ``minNetMILP``.** Its binary formulation forces a kept reaction to carry
at least ``tol`` and a dropped one to carry exactly zero, so the kept flux
vector is exactly mass balanced and the question does not arise. The cost is
a genome-scale MILP.

.. code-block:: python

   extracted = spectra_me(model, core, tol=1e-4, problem_type="minNetMILP",
                          time_limit=900, seed=0)

**Or prune afterwards**, accepting that core reactions may be lost:

.. code-block:: python

   from spectra_cobra import consistent_reaction_ids

   pruned = extracted.copy()
   while True:
       ids = {r.id for r in pruned.reactions}
       ok, _ = consistent_reaction_ids(pruned, tol=1e-4, seed=0)
       if ids == ok:
           break
       pruned.remove_reactions(sorted(ids - ok), remove_orphans=True)

This converges in a couple of rounds and does give a consistent model, but it
removed between two and five core reactions in the trials above — so check
what survived before relying on it:

.. code-block:: python

   survived = {r.id for r in pruned.reactions}
   print(f"core kept after pruning: {sum(c in survived for c in core)}/{len(core)}")

Raising ``tol`` does not help, and lowering it makes the dust problem worse.

Occasional direction failures
-----------------------------

One of the twelve trials raised:

.. code-block:: text

   SpectraError: 2 core reaction(s) could not be given a flux direction ...
   Their accumulated flux stays below 1e-05, which can happen when the convex
   combination of two iterations cancels out; retry with a different seed.

The direction phase blends successive flux vectors as a convex combination, so
two iterations can cancel on a given reaction and leave it at zero. It is a
property of the algorithm, not of your core set, and a different ``seed``
fixes it. Worth wrapping if you are running many core sets unattended:

.. code-block:: python

   from spectra_cobra import SpectraError

   for attempt in range(5):
       try:
           extracted = spectra_me(model, core, tol=1e-4, seed=attempt)
           break
       except SpectraError:
           continue
   else:
       raise RuntimeError("no seed worked for this core set")

Performance notes
-----------------

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Step
     - Recon3D (11303 rxns)
     - Textbook (95 rxns)
   * - ``spectra_cc``
     - ~15 s, 8 LPs
     - ~0.1 s, 8 LPs
   * - ``spectra_me`` (``minNetLP``)
     - ~50–75 s
     - ~0.05 s
   * - ``spectra_me`` (``minNetMILP``)
     - minutes; set ``time_limit``
     - ~0.8 s

Most of the ``minNetLP`` time is spent constructing the LP — one
absolute-value variable and two constraints per candidate reaction, so around
11000 variables and 22000 constraints — not solving it. If you are extracting
many models from one universal model, that construction cost is paid per call.

Reproducing this
----------------

The scripts behind this page are not shipped with the package, since Recon3D
is not redistributed here. Get the model from the
`MATLAB repository <https://github.com/bisect-group/spectra>`_ at
``Recon3D+/consRecon3DGeneSymbol.mat``. Note that cobrapy's
``load_matlab_model`` rejects it, because the GeneSymbol variant has duplicate
gene symbols; load the stoichiometry with :mod:`scipy.io` and build the model
without gene annotation.

Next
----

* :doc:`formulations` — the trade-offs between the formulations
* :doc:`../matlab_differences` — what else differs from MATLAB
