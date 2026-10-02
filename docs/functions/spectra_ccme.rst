spectra_ccme
============

Checks consistency and extracts a context-specific model in a single pass.
Ports ``spectraCCME.m``.

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

Reference
---------

.. autofunction:: spectra_cobra.spectra_ccme

See also
--------

* :doc:`spectra_me` — when the model is already consistent
* :doc:`../tutorials/extraction` — worked examples
