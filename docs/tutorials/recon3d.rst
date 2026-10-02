Genome scale: Recon3D
=====================

Everything so far used toy models and the 95-reaction textbook model. This
page runs the same workflow on Recon3D, the human reconstruction SPECTRA's own
experiments use, and reports what actually happens — including where the
default formulation falls short.

The model used here is the consistent Recon3D from the MATLAB repository:
11303 reactions over 6388 metabolites.

Getting a fully consistent model
--------------------------------

.. code-block:: python

   from spectra_cobra import spectra_cc

   model.solver = "gurobi"
   consistent = spectra_cc(model, tol=1e-4, detection_cutoff=1e-7, seed=0)
   print(f"{len(model.reactions)} -> {len(consistent.reactions)} reactions")

On Recon3D restricted to a defined medium — which is where context-specific
modelling actually starts, most exchanges closed — that gives:

.. code-block:: text

   11303 -> 8406 reactions

2897 reactions removed, in 8 LPs and about 30 seconds.

.. note::

   The files in the MATLAB repository's ``Recon3D+/`` directory are *already*
   consistency-checked: both ``UpdatedRecon3D.mat`` and
   ``consRecon3DGeneSymbol.mat`` hold 11303 reactions with nothing blocked, as
   does the ``Reconmodel.mat`` used in the PCOS study (10600 reactions, 0
   blocked). Running the check on one of those tells you only that the file is
   what it claims. Close the medium first, or start from an unchecked
   reconstruction, if you want to see the check do work.

Proving it is consistent
------------------------

A fully consistent model is a **fixed point**: running the check again must
remove nothing. Use a different seed, so you are not just repeating the same
randomised tie-breaking:

.. code-block:: python

   from spectra_cobra import consistent_reaction_ids

   again, n_lps = consistent_reaction_ids(consistent, tol=1e-4,
                                          detection_cutoff=1e-7, seed=1)
   assert len(again) == len(consistent.reactions), "not actually consistent"

.. code-block:: text

   8406/8406 consistent in 8 LPs

And confirm independently with cobrapy, which uses FVA rather than this
algorithm:

.. code-block:: python

   from cobra.flux_analysis import find_blocked_reactions

   assert find_blocked_reactions(consistent, zero_cutoff=1e-7) == []

Both hold. The same cycle on iJO1366, which ships genuinely inconsistent,
removes 878 of 2583 reactions and the remaining 1705 are a fixed point with
zero blocked by FVA.

Accuracy against fastcc
-----------------------

On the defined-medium Recon3D the two methods disagree on six reactions.
Taking FVA on the full model as the arbiter:

.. list-table::
   :header-rows: 1
   :widths: 30 20 25 25

   * - Reaction
     - FVA \|v\|max
     - Truly
     - Which was right
   * - ``EX_M00234[e]``
     - 0
     - blocked
     - ``spectra_cc``
   * - ``EX_M01807[e]``
     - 0
     - blocked
     - ``spectra_cc``
   * - ``HMR_0002``
     - 0
     - blocked
     - ``spectra_cc``
   * - ``HMR_3996``
     - 1000
     - consistent
     - ``spectra_cc``
   * - ``r1391``
     - 1000
     - consistent
     - ``spectra_cc``
   * - ``r1392``
     - 1000
     - consistent
     - ``spectra_cc``

``spectra_cc`` is right on all six; ``fastcc`` produced three false positives
and three false negatives. On the PCOS ``Reconmodel.mat`` the pattern repeats
on a smaller scale: ``spectra_cc`` and FVA agree that all 10600 reactions are
consistent, while ``fastcc`` reports 3 of them blocked.

Reproducibility at genome scale
-------------------------------

Part of that disagreement is that ``fastcc`` does not give the same answer
twice. Its singleton phase picks an arbitrary element from a Python set, and
set iteration order depends on ``PYTHONHASHSEED``, so the result shifts from
run to run. ``spectra_cc`` given a ``seed`` does not:

.. list-table::
   :header-rows: 1
   :widths: 34 33 33

   * - Run
     - ``spectra_cc`` (``seed=0``)
     - ``fastcc``
   * - ``PYTHONHASHSEED=1``
     - 8406
     - 8408
   * - ``PYTHONHASHSEED=2``
     - 8406
     - 8405
   * - ``PYTHONHASHSEED=3``
     - 8406
     - 8403
   * - two further runs
     - 8406, 8406
     - 8410, 8406

Across five runs on the same model ``fastcc`` spans 8403 to 8410, a range of
seven reactions, while ``spectra_cc`` returns 8406 every time.

This is not a claim that ``spectra_cc`` is always more accurate — six
reactions out of 11303 is a 0.05% disagreement, and both are far closer to
each other than to being wrong. But it does mean the port is not merely
reproducing ``fastcc``, that it is the reproducible of the two, and that at
genome scale a consistency result is worth confirming with FVA whichever
method produced it.

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

   The cause is the inclusion rule, which comes from MATLAB: a reaction is
   kept if its flux in the extraction LP exceeds ``tol * 1e-7``, floored here
   at ``model.tolerance``, so around ``1e-7``. A handful of reactions end up
   carrying flux *just under* that cutoff while doing load-bearing balancing
   work. Dropping them leaves a mass-balance residual of the same order — in
   the trial examined, up to ``3.5e-07`` on ``coa[c]``, on 3 of 972
   metabolites — which is above the solver's feasibility tolerance, so the
   kept flux vector is not actually feasible in the extracted model. For the
   chains that relied on those sub-cutoff reactions the only feasible flux is
   then exactly zero: all 25 sampled blocked reactions had an attainable flux
   of **0**, despite carrying around ``1e-4`` in the LP. They were mostly
   exchange and transport pairs, such as ``EX_dxtrn[e]`` with ``DXTRNt``.

   Note what this is *not*. The total flux discarded is tiny — 1.6e-06 across
   all 10376 dropped reactions — so this is not an accumulation of noise. It
   is a few specific reactions below the cutoff whose absence breaks a chain
   outright.

   **Always check the result** rather than assuming it:

   .. code-block:: python

      from spectra_cobra import blocked_reaction_ids

      blocked = set(blocked_reaction_ids(extracted, tol=1e-4))
      assert not blocked, f"{len(blocked)} blocked, {len(blocked & set(core))} core"

Getting a consistent model
--------------------------

Two options, with different costs.

**Use ``minNetMILP``.** A binary at zero forces its reaction's flux to
*exactly* zero, so nothing sub-cutoff is discarded and the kept flux vector is
exactly mass balanced. That is an argument from the formulation rather than a
measurement: at a core size of 50 the genome-scale MILP did not finish within
15 minutes on Gurobi, so it is untested at this scale, and the cost is real.
Set ``time_limit`` and expect a feasible-but-not-proven-optimal answer.

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
