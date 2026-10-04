spectra_me
==========

Extracts a context-specific model around a set of core reactions.

The input model must be flux consistent. Run :doc:`spectra_cc` first, or use
:doc:`spectra_ccme`, which folds the check in.

Usage
-----

.. code-block:: python

   from spectra_cobra import spectra_me

   extracted = spectra_me(model, core_reactions=["PGI", "PFK"], tol=1e-4)

How it works
------------

Two phases.

**Phase 1, settling directions.** The forward and reverse LPs alternate over
the core reactions, accumulating a flux vector whose signs say which direction
each core reaction should run in. Successive solutions are blended as a convex
combination, :math:`c \cdot v_{old} + (1 - c) \cdot v_{new}` with :math:`c`
drawn from [0.45, 0.55], so a direction preference carries across iterations.

**Phase 2, network inference.** One optimisation problem picks the rest of the
network subject to those directions. Which problem is up to you; see
:doc:`formulations`.

Arguments worth knowing
-----------------------

``core_reactions``
   Accepts reaction identifiers, :class:`cobra.Reaction` objects, or integer
   indices into ``model.reactions``.

``weights``
   A dict keyed by reaction identifier; reactions left out default to 1.0.
   What a weight *means* depends on ``problem_type`` — see
   :doc:`formulations`.

``problem_type``
   One of ``"minNetLP"`` (default), ``"minNetMILP"``, ``"tradeOff"`` or
   ``"growthOptim"``.

``n_solutions`` and ``alt_solution_method``
   For several alternative models; see
   :doc:`../tutorials/alternative_solutions`.

``seed``
   Pins the randomised objective coefficients so a run is reproducible.

Errors
------

:class:`~spectra_cobra.SpectraInfeasibleCoreError`
   No flux distribution lets every core reaction carry ``tol`` at once. Either
   the core set is mutually inconsistent, or some of its reactions are blocked
   in the model. :doc:`spectra_ccme` handles the latter by dropping them.

Reference
---------

.. autofunction:: spectra_cobra.spectra_me

See also
--------

* :doc:`../tutorials/extraction` — worked examples
* :doc:`../tutorials/recon3d` — a genome-scale walkthrough
* :doc:`formulations` — choosing what to optimise
