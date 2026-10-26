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

The consistency loop of :doc:`spectra_cc` runs first, but the flux it finds is
kept rather than discarded, so the core directions are read off it instead of
being recomputed. That saves a pass over the model compared with running the
two functions in sequence. The network inference step is then the same as in
:doc:`spectra_me`.

Arguments worth knowing
-----------------------

Everything :doc:`spectra_me` takes is accepted here too, with four
exceptions noted at the end.

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
   ``"growthOptim"``.

``n_solutions`` and ``alt_solution_method``
   For several alternative models; see
   :doc:`../tutorials/alternative_solutions`. ``alt_solution_method`` is
   ``"coreDirection"`` by default, with ``"pathwayExclusion"`` available for
   the mixed-integer formulations.

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

Not accepted here
-----------------

``indicator_reactions``, ``indicator_groups``, ``previous_solutions`` and
``return_solutions`` belong to :doc:`spectra_me` alone. If you need them on a
model that may be inconsistent, run :doc:`spectra_cc` first and then
:doc:`spectra_me`, which is what this function does internally anyway.

Reference
---------

.. autofunction:: spectra_cobra.spectra_ccme

See also
--------

* :doc:`spectra_me` — when the model is already consistent
* :doc:`../tutorials/extraction` — worked examples
