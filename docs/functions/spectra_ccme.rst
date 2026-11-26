spectra_ccme
============

Checks consistency and extracts a context-specific model in a single pass.

Use this instead of :doc:`spectra_cc` followed by :doc:`spectra_me` when the
universal model is not known to be flux consistent. Blocked core reactions
are dropped and reported rather than causing a failure.

Usage
-----

.. code-block:: python

   from spectra_cobra import spectra_ccme

   extracted, blocked_core = spectra_ccme(model, core_reactions=["PGI", "PFK"],
                                          tol=1e-4)
   if blocked_core:
       print(f"dropped, because blocked: {blocked_core}")

Note the two return values: the model, and the list of core reaction
identifiers that turned out to be blocked and so are absent from it.

How it works
------------

Stage 1a prunes structural dead ends in 0 LPs, then runs a core-scoped
pinned ``forward_cc``/``reverse`` loop (sharing one random weight draw
:math:`w_i \sim U(1, 1.1)` per round, with zero convex combinations) together
with an unpinned verification guard so mutually conflicting finite-bound core
reactions emit a warning rather than being misreported in ``blocked_core``.
Stage 1b refines two-way active core signs with a 1-LP Anchor when applicable,
and Stage 2 dispatches to the chosen formulation with consistency certification,
returning a copy with original cobrapy reaction orientations and bounds
preserved.

Arguments worth knowing
-----------------------

Everything :doc:`spectra_me` takes (except ``previous_solutions``) is accepted
here too.

``core_reactions``
   Accepts reaction identifiers, :class:`cobra.Reaction` objects, or integer
   indices into ``model.reactions``. Any that turn out to be blocked come
   back in the second return value rather than raising.

``consistency_type``
   ``"stoichiometry"`` (default) for a steady state, or ``"topology"`` for
   the accumulation condition. This governs both halves of the call: the
   consistency check and the extraction. See :doc:`../tutorials/topology`.

``weights``
   A dict keyed by reaction identifier; reactions left out default to 1.0.
   What a weight *means* depends on ``problem_type`` — see
   :doc:`formulations`.

``problem_type``
   One of ``"minNetLP"`` (default), ``"minNetMILP"``, ``"tradeOff"`` or
   ``"optimBiomass"``.

``n_solutions`` and ``alt_solution_method``
   For several alternative models; see
   :doc:`../tutorials/alternative_solutions`. ``alt_solution_method`` defaults
   to ``"coreDirection"`` for LP formulations and ``"pathwayExclusion"`` for
   MILP formulations.

``time_limit``
   Seconds to spend on a mixed-integer solve (default 7200). A solution
   found before the cap is returned even if not proven optimal.

``remove_genes``
   Whether to drop genes left without a reaction (default True). Orphaned
   metabolites are always dropped.

``seed``
   Pins the randomised objective coefficients so a run is reproducible.

``inclusion_cutoff``
   The absolute flux at which a reaction counts as part of the extracted
   model (default ``tol * 1e-7``). The default sits far below the solver's
   own tolerance deliberately: raising it drops reactions that carry almost
   no flux but are load-bearing for a mass balance, which leaves core
   reactions in the result yet unable to carry flux.

   For the LP formulations this is the whole membership rule. The
   mixed-integer formulations take the **union** of their binaries and this
   rule, so a reaction can be in the result with its binary switched off,
   having got there on flux alone. It is not ignored there — see
   :doc:`../tutorials/recon3d`, where that union is what made ``minNetMILP``
   appear to return larger models than ``minNetLP``.

``feasibility_tolerance``, ``optimality_tolerance``, ``integrality_tolerance``
   The solver tolerances both halves of the pass run at, each defaulting to
   ``1e-9``, the lowest most solvers accept. The consistency check decides
   whether a reaction reached ``tol``, and the extraction reads its membership
   off a solution, so a loose feasibility tolerance can let a reaction pass
   the check on residual rather than real flux and then sit in the result
   unable to carry any. Pass ``None`` for any of them to leave the solver's
   own in place. ``integrality_tolerance`` reaches only the mixed-integer
   ``problem_type`` values; neither the consistency check nor the linear
   formulations have binaries.

``indicator_reactions``, ``indicator_groups``, and ``return_solutions``
   Supported directly on :func:`~spectra_cobra.spectra_ccme` for MILP group
   selection and solution inspection. To pass explicit ``previous_solutions``,
   run :doc:`spectra_cc` followed by :doc:`spectra_me`.

Reference
---------

.. autofunction:: spectra_cobra.spectra_ccme

See also
--------

* :doc:`spectra_me` — when the model is already consistent
* :doc:`../tutorials/extraction` — worked examples
