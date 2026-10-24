Genome scale: Recon3D
=====================

Everything so far used toy models and the 95-reaction textbook model. This
page runs the same workflow on Recon3D, the human reconstruction SPECTRA's own
experiments use, and reports what actually happens — including where the
default formulation falls short.

Get the Recon3D model from the COBRA Toolbox's model collection:

   https://github.com/opencobra/COBRA.models/blob/master/mat/Recon3DModel_301.mat

That is ``Recon3DModel_301.mat`` — **10600 reactions** over 5835 metabolites
and 2248 genes, with ``biomass_maintenance`` as its objective. It is the same
model the COBRA Toolbox ships for its own tests, so it needs no conversion
and no registration. Everything below was measured against it.

.. code-block:: python

   from cobra.io import load_matlab_model

   model = load_matlab_model("Recon3DModel_301.mat")
   model.solver = "gurobi"
   print(f"{len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

.. code-block:: text

   10600 reactions, 5835 metabolites

.. note::

   Other Recon3D distributions may not load as cleanly. Some MATLAB
   structures carry a duplicated gene — ``AOC2`` is the usual culprit — which
   older cobrapy versions refuse. If you hit that, convert with the COBRA
   Toolbox or take Recon3D from
   `BiGG <http://bigg.ucsd.edu/models/Recon3D>`_. The file linked above
   loads directly.

Watching the check do work
--------------------------

Recon3D as distributed is *already* consistency-checked — the 10600-reaction
file above has nothing blocked, as the measurement below shows. Running the
check on a model like that tells you only that the file is what it claims.

To see it remove anything you need a reconstruction that ships inconsistent.
iJO1366 does, and it comes with cobrapy, so this runs anywhere:

.. code-block:: python

   from cobra.io import load_model
   from spectra_cobra import spectra_cc

   model = load_model("iJO1366")
   model.solver = "gurobi"
   consistent = spectra_cc(model, tol=1e-4, detection_cutoff=1e-7, seed=0)
   print(f"{len(model.reactions)} -> {len(consistent.reactions)} reactions")

.. code-block:: text

   2583 -> 1705 reactions

878 reactions removed, in about 1.5 seconds.

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

   1705/1705 consistent in 12 LPs

And confirm independently with cobrapy, which uses FVA rather than this
algorithm:

Both hold: the 1705 kept reactions are a fixed point.

.. note::

   Use :func:`~spectra_cobra.spectra_cc` for this rather than cobrapy's
   ``find_blocked_reactions``. That function is FVA underneath — two LPs per
   reaction — which at genome scale takes long enough to be impractical, and
   cobrapy's ``fastcc`` is not a reliable substitute. Running the consistency
   check again and seeing nothing removed answers the same question in one
   pass: on the 10600-reaction Recon3D it takes about 6 seconds.

Recon3D itself is already consistent, so the check removes nothing:

.. code-block:: python

   consistent_recon3d = spectra_cc(model, tol=1e-4, detection_cutoff=1e-7,
                                   seed=0)
   print(f"{len(model.reactions)} -> {len(consistent_recon3d.reactions)}")

   # and it is a fixed point under a different seed
   again = spectra_cc(consistent_recon3d, tol=1e-4, detection_cutoff=1e-7,
                      seed=1)
   print(f"re-check: {len(again.reactions)}/{len(consistent_recon3d.reactions)}")

.. code-block:: text

   10600 -> 10600
   re-check: 10600/10600

Which is the point made above: running the check on a model that is already
consistent tells you only that the file is what it claims.

Extraction with a random core set
---------------------------------

A realistic test: draw a random core set, extract, and check the two things
the method promises — that every core reaction is in the result, and that the
result is flux consistent.

Draw the core from the *consistent* model and extract from it.
:func:`~spectra_cobra.spectra_me` requires a consistent input, and on Recon3D
that costs nothing, since the model is already consistent.

.. code-block:: python

   import numpy as np
   from spectra_cobra import spectra_cc, spectra_me

   pool = sorted(r.id for r in consistent_recon3d.reactions)

   for size in (10, 50, 200):
       for trial in range(3):
           rng = np.random.default_rng(trial)
           core = sorted(rng.choice(pool, size=size, replace=False).tolist())

           extracted = spectra_me(consistent_recon3d, core, tol=1e-4,
                                  problem_type="minNetLP", seed=0)
           kept = {rxn.id for rxn in extracted.reactions}

           # Verify with the consistency check, not with FVA: a consistent
           # result loses nothing when the check is run on it again.
           checked = spectra_cc(extracted, tol=1e-4, detection_cutoff=1e-7,
                                seed=2)
           blocked = len(extracted.reactions) - len(checked.reactions)

           print(f"core {size:>3} trial {trial}: {len(kept):>4} reactions, "
                 f"core kept {sum(c in kept for c in core)}/{size}, "
                 f"blocked {blocked}")

.. code-block:: text

   core  10 trial 0:  136 reactions, core kept 10/10, blocked 0
   core  10 trial 1:  221 reactions, core kept 10/10, blocked 0
   core  10 trial 2:  109 reactions, core kept 10/10, blocked 0
   core  50 trial 0:  460 reactions, core kept 50/50, blocked 0
   core  50 trial 1:  508 reactions, core kept 50/50, blocked 0
   core  50 trial 2:  444 reactions, core kept 50/50, blocked 0
   core 200 trial 0: 1286 reactions, core kept 200/200, blocked 0
   core 200 trial 1: 1246 reactions, core kept 200/200, blocked 0
   core 200 trial 2: 1298 reactions, core kept 200/200, blocked 0

.. list-table::
   :header-rows: 1
   :widths: 14 14 26 24 22

   * - core
     - trials
     - extracted
     - core kept
     - blocked
   * - 10
     - 3
     - 109–221
     - all
     - 0
   * - 50
     - 3
     - 444–508
     - all
     - 0
   * - 200
     - 3
     - 1286–1298
     - all
     - 0

Both promises hold in every trial, and each one takes five to seven seconds
including the verification. The spread within a core size is wide — 109
against 221 reactions for ten core reactions — because what a core set costs
depends on where in the network it falls, not on how large it is.

.. note::

   Earlier versions of this page reported blocked reactions in most trials,
   including blocked *core* reactions. That came from two things, both since
   fixed: the core was being drawn from, and extracted out of, the full model
   rather than the consistent one, and ``inclusion_cutoff`` used to be
   floored at ``model.tolerance``. The next section explains the second, which
   is worth understanding even though the default now handles it, because the
   same failure returns if you raise the cutoff yourself.

Two settings that decide whether it works
-----------------------------------------

The blocked core reactions this page used to report came from two numerical
settings, not from the algorithm. Both are handled by the defaults now, and
the trials above show none — but raising either one yourself brings the
failure back, so it is worth knowing what they do.

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

Verify with :func:`~spectra_cobra.spectra_cc`, by running the consistency
check on the extracted model and seeing that it removes nothing:

.. code-block:: python

   checked = spectra_cc(extracted, tol=1e-4, detection_cutoff=1e-7, seed=2)
   blocked = len(extracted.reactions) - len(checked.reactions)
   assert blocked == 0, f"{blocked} reactions in the result cannot carry flux"

:func:`~spectra_cobra.check_extraction` wraps this up with the core-set check
as well. It still defaults to ``method="fva"``; pass ``method="spectra"`` at
genome scale.

.. code-block:: python

   from spectra_cobra import check_extraction

   report = check_extraction(extracted, core, tol=1e-4, method="spectra")
   assert not report.blocked_core, report.summary()

.. note::

   **Avoid cobrapy's ``find_blocked_reactions`` at this size.** It is FVA
   underneath — two LPs for every reaction — where the consistency check
   converges in a handful of passes over the whole model. On the
   1286-reaction extraction above the two agree exactly, both reporting no
   blocked reactions, but FVA takes 4 seconds against 0.4, and the gap widens
   with the model. cobrapy's ``fastcc`` is not a dependable substitute
   either.

   An earlier version of this page recommended the opposite, because the
   LP-driven check had reported no blocked core reactions where FVA found
   four. That was the ``inclusion_cutoff`` behaviour described above, which
   left numerically marginal reactions in the result for the check to be
   fooled by. With the current default the two agree, and the faster one is
   the one to use.

   If FVA does raise a bare ``RuntimeError``, that is its process pool, not
   your model: call it from inside an ``if __name__ == "__main__":`` guard,
   or pass ``processes=1``.

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

Each row below gave Gurobi ten minutes on an 11303-reaction consistent
Recon3D variant — a different file from the 10600-reaction one used above,
so read the rows against each other rather than against the trials — against
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

Recon3D is not redistributed with this package. Get
``Recon3DModel_301.mat`` from the COBRA Toolbox's model collection:

   https://github.com/opencobra/COBRA.models/blob/master/mat/Recon3DModel_301.mat

Everything on this page except the ``minNetMILP`` comparison was measured
against that file, with Gurobi, and the code shown is the code that produced
the numbers.

References
----------

* Brunk, E., Sahoo, S., Zielinski, D. C., et al. (2018). Recon3D enables a
  three-dimensional view of gene variation in human metabolism. *Nature
  Biotechnology*, 36(3), 272–281. https://doi.org/10.1038/nbt.4072

Next
----

* :doc:`formulations` — the trade-offs between the formulations
