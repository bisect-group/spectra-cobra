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

**Phase 1, discovering core directions.** Core reactions are first classified
in 0 LPs by bound reachability at ``tol``: empty core sets and all-one-way
core sets skip 100% of Stage-1 LPs, and bound-unreachable core reactions raise
:class:`~spectra_cobra.SpectraInfeasibleCoreError` immediately. When two-way
reversible core reactions exist alongside fixed flux, a single Pass-0
``minNetLP`` Anchor LP steers directions along low-cost (or non-penalised
``tradeOff``) pathways, and any remaining unsettled two-way core reactions are
settled by iterative pinned ``forward``/``reverse`` LPs sharing one random
weight draw :math:`w_i \sim U(1, 1.1)` per round (with zero convex
combinations). Reverse-only irreversible reactions (``upper_bound <= 0``) are
oriented internally during optimisation and returned with their original
cobrapy stoichiometry and bounds intact.

**Phase 2, network inference.** One optimisation problem picks the rest of the
network subject to those directions (imposing :math:`|v_i| \ge \text{tol}`
simultaneously across all core reactions along their settled directions). Which
problem is up to you; see :doc:`formulations`.

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
   ``"optimBiomass"``.

``n_solutions`` and ``alt_solution_method``
   For several alternative models; see
   :doc:`../tutorials/alternative_solutions`. When ``alt_solution_method`` is
   ``None`` (the default), selects ``"coreDirection"`` for LP formulations and
   ``"pathwayExclusion"`` for MILP formulations.

``seed``
   Pins the randomised objective coefficients so a run is reproducible.

``indicator_reactions`` and ``indicator_groups``
   Change what the mixed-integer binaries range over: one per named
   reaction, or one shared between the reactions of each group. See
   :doc:`formulations`.

``return_solutions``
   Also hand back the :class:`~spectra_cobra.MilpSolution` behind each model.
   The solution says which indicators the solver actually switched on, which
   the model cannot: a reaction can be in the model because it carries a
   trace of flux while its binary is off. When the binaries mean something in
   themselves — an organism's presence, a group's membership — that
   distinction is the answer rather than an implementation detail.

   .. code-block:: python

      model, solution = spectra_me(..., return_solutions=True)
      solution.selected      # the indicators switched on
      solution.included      # those, plus whatever carries flux regardless

``feasibility_tolerance``, ``optimality_tolerance``, ``integrality_tolerance``
   The solver tolerances the extraction runs at, each defaulting to ``1e-9``,
   the lowest most solvers accept. These are not cosmetic: the extraction
   reads its answer off a solution rather than off the objective, so a loose
   feasibility tolerance leaves core reactions present but unable to carry
   flux, and a loose integrality tolerance keeps reactions whose binary the
   solver switched off. Both failures are silent. Pass ``None`` for any of
   them to leave the solver's own in place. See
   :doc:`../tutorials/recon3d` for what each one costs when it is wrong.

   ``integrality_tolerance`` reaches only the mixed-integer ``problem_type``
   values; ``minNetLP`` and ``optimBiomass`` have no binaries for it to act on.

Thermodynamic cone extraction (``spectra_me2``)
-----------------------------------------------

:func:`~spectra_cobra.spectra_me2` extracts a subnetwork on the pure
thermodynamic cone — replacing finite flux bounds with ``(-inf, inf)``,
``(0, inf)``, or ``(-inf, 0)`` according to reaction reversibility — while
preserving original cobrapy reaction orientations.

Errors
------

:class:`~spectra_cobra.SpectraInfeasibleCoreError`
   No flux distribution lets every core reaction carry ``tol`` at once. Either
   the core set is mutually inconsistent, or some of its reactions are blocked
   in the model. :doc:`spectra_ccme` handles the latter by dropping them.

Reference
---------

.. autofunction:: spectra_cobra.spectra_me

.. autofunction:: spectra_cobra.spectra_me2

See also
--------

* :doc:`../tutorials/extraction` — worked examples
* :doc:`../tutorials/recon3d` — a genome-scale walkthrough
* :doc:`formulations` — choosing what to optimise
