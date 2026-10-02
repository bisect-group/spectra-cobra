API reference
=============

A flat index of everything the package exports. The pages under *Functions*
cover the same material with explanation and worked examples; this is for
looking up a signature.

Consistency
-----------

.. autofunction:: spectra_cobra.spectra_cc
   :no-index:
.. autofunction:: spectra_cobra.consistent_reaction_ids
   :no-index:
.. autofunction:: spectra_cobra.blocked_reaction_ids
   :no-index:

Extraction
----------

.. autofunction:: spectra_cobra.spectra_me
   :no-index:
.. autofunction:: spectra_cobra.spectra_ccme
   :no-index:

Formulations
------------

.. autofunction:: spectra_cobra.min_net_lp
   :no-index:
.. autofunction:: spectra_cobra.min_net_milp
   :no-index:
.. autofunction:: spectra_cobra.trade_off
   :no-index:
.. autofunction:: spectra_cobra.growth_optim
   :no-index:

Flux
----

.. autofunction:: spectra_cobra.flux_reducer
   :no-index:

Verification
------------

.. autofunction:: spectra_cobra.check_extraction
   :no-index:
.. autoclass:: spectra_cobra.ExtractionReport
   :no-index:
   :no-members:

Exceptions
----------

.. autoexception:: spectra_cobra.SpectraError
   :no-index:
.. autoexception:: spectra_cobra.SpectraSolverError
   :no-index:
.. autoexception:: spectra_cobra.SpectraInfeasibleCoreError
   :no-index:
