Genome scale: Recon3D
=====================

Everything so far used toy models and the 95-reaction textbook model. This
page runs the same workflow on Recon3D, the human reconstruction SPECTRA's own
experiments use, and reports what actually happens — including where the
default formulation falls short.

The model used here is the consistent Recon3D distributed with SPECTRA:
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

   The files in the SPECTRA repository's ``Recon3D+/`` directory are *already*
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

   The cause was the inclusion rule. A reaction is kept if its flux in the
   extraction LP exceeds ``inclusion_cutoff``, and that cutoff used to be
   floored at ``model.tolerance``, so around ``1e-7``. A handful of reactions
   carry flux *just under* that while doing load-bearing balancing work.
   Dropping them leaves a mass-balance residual of the same order — up to
   ``3.5e-07`` on ``coa[c]``, on 3 of 972 metabolites — above the solver's
   feasibility tolerance, so the kept flux vector is not actually feasible in
   the extracted model. The chains that relied on them then collapse: all 25
   sampled blocked reactions had an attainable flux of **0** despite carrying
   around ``1e-4`` in the LP, and they were mostly exchange and transport
   pairs such as ``EX_dxtrn[e]`` with ``DXTRNt``.

   Note what this was *not*: the total flux discarded is tiny, 1.6e-06 across
   all 10376 dropped reactions, so it was never an accumulation of noise. It
   was a few specific sub-cutoff reactions whose absence breaks a chain
   outright.

Two settings that decide whether it works
-----------------------------------------

The blocked core reactions above came from two numerical settings, not from
the algorithm. Fixing both removes the problem entirely.

**1. Tighten the solver.** This is the big one, and the least obvious. The
extraction reads its answer off an LP solution whose mass balance holds only
to the solver's feasibility tolerance. At cobrapy's default of ``1e-7``, that
residual is enough to break chains carrying flux of order ``tol = 1e-4``. Set
it as low as the solver allows — ``1e-9`` is the floor on Gurobi and CPLEX:

.. code-block:: python

   model.tolerance = 1e-9

**2. Do not floor the inclusion cutoff.** A reaction joins the extracted
model when its flux clears ``inclusion_cutoff``, which defaults to
``tol * 1e-7`` — far below the solver tolerance, deliberately. Raising it to
something that looks more sensible discards reactions that carry almost no
flux but are load-bearing for a mass balance, and the chains that relied on
them collapse.

Measured over the same twelve Recon3D trials, with every core reaction's
attainable flux checked by FVA:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Configuration
     - Dead core reactions
     - Trials affected
   * - cutoff floored at ``model.tolerance``, solver ``1e-7``
     - 16
     - 5 of 12, plus 1 hard failure
   * - cutoff ``tol * 1e-7``, solver ``1e-7``
     - 4
     - 2 of 12
   * - cutoff ``tol * 1e-7``, solver ``1e-9``
     - **0**
     - **none**

With both in place every one of the twelve trials returns a model containing
all of its core reactions, every one of them able to carry flux, with no
failures — and slightly faster, at 46–50 s rather than 50–61 s.

:func:`~spectra_cobra.spectra_me` warns when ``model.tolerance`` is within a
factor of 1e4 of ``tol``, since this is not a setting most callers would
think to check.

.. code-block:: text

   model.tolerance is 1e-07 against tol=0.0001, a ratio of only 1000. The
   extraction reads its answer off an LP whose mass balance holds to the
   solver's tolerance, so a loose one can leave core reactions present in the
   result but unable to carry flux. Set model.tolerance=1e-9, the lowest most
   solvers accept, before extracting.

Verifying the result
--------------------

.. warning::

   **Verify with FVA, not with this package's own check.** On these
   extractions ``consistent_reaction_ids`` reported no blocked core reactions
   where FVA found four. The extracted models hold numerically marginal
   reactions by design, and the LP-driven check can call a dead one
   consistent. :func:`~spectra_cobra.check_extraction` therefore defaults to
   ``method="fva"``:

   .. code-block:: python

      from spectra_cobra import check_extraction

      report = check_extraction(extracted, core, tol=1e-4)   # FVA by default
      assert not report.blocked_core, report.summary()

   ``method="spectra"`` is available and much faster, but treat its verdict as
   a screen rather than a guarantee.

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
     - does not finish; set ``time_limit``
     - ~0.8 s

Most of the ``minNetLP`` time is spent constructing the LP — one
absolute-value variable and two constraints per candidate reaction, so around
11000 variables and 22000 constraints — not solving it. If you are extracting
many models from one universal model, that construction cost is paid per call.

Is minNetMILP worth it at this scale?
-------------------------------------

``minNetLP`` minimises the L1 norm of the flux as a proxy for the number of
reactions; ``minNetMILP`` counts them exactly, with one binary per candidate.
At genome scale that is an 11000-binary problem, and Gurobi does not come
close to solving it. It is still worth running, because the proxy is loose
enough that even a poor incumbent beats the LP by a wide margin.

Each row below gave Gurobi ten minutes on the 11303-reaction model, against
``minNetLP`` on the same core set:

.. list-table::
   :header-rows: 1
   :widths: 10 18 20 14 14 24

   * - Core
     - ``minNetLP``
     - ``minNetMILP``
     - Smaller by
     - MIP gap
     - Shared / LP only / MILP only
   * - 10
     - 170
     - 100
     - 41%
     - 65.6%
     - 78 / 92 / 22
   * - 50
     - 662
     - 439
     - 34%
     - 53.5%
     - 373 / 289 / 66
   * - 200
     - 1791
     - 1551
     - 13%
     - 36.9%
     - 1388 / 403 / 163

Every model kept its whole core set, and no core reaction was blocked in any
of them.

Two things are worth reading off this. The advantage shrinks as the core
grows, because a larger core forces more of the network and leaves less to
choose; by a core of 200 the LP is within 13%. And the mixed-integer model is
not a pruned version of the linear one — it includes reactions the LP left
out, so the two pick genuinely different routes rather than one being a
subset of the other.

The solutions are far from proven: at a core of 10 the bound was 31 against
an incumbent of 90, so a much smaller model may well exist. Treat
``minNetMILP`` at this scale as "a better answer than the LP for ten minutes
of effort", not as the minimum network.

.. note::

   A solve stopped at ``time_limit`` returns its incumbent and logs a
   warning. Give it more time if the gap matters; the limit binds to the
   second, so budget it directly.

Reproducing this
----------------

The scripts behind this page are not shipped with the package, since Recon3D
is not redistributed here. Get the model from the
`SPECTRA repository <https://github.com/bisect-group/spectra>`_ at
``Recon3D+/consRecon3DGeneSymbol.mat``. Note that cobrapy's
``load_matlab_model`` rejects it, because the GeneSymbol variant has duplicate
gene symbols; load the stoichiometry with :mod:`scipy.io` and build the model
without gene annotation.

Next
----

* :doc:`formulations` — the trade-offs between the formulations
