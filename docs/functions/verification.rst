check_extraction
================

Checks an extracted model against what the extraction promised: that every
core reaction is present, and that the result is flux consistent.

Usage
-----

.. code-block:: python

   from spectra_cobra import spectra_me, check_extraction

   extracted = spectra_me(model, core, tol=1e-4)
   report = check_extraction(extracted, core, tol=1e-4)

   if not report.is_valid:
       print(report.summary())

Pass the same ``tol`` and ``consistency_type`` the extraction used, or the
report will not describe the model you actually have.

Why this exists
---------------

The first promise always holds; the second does not.

``minNetLP`` and ``growthOptim`` read their answer off a flux vector, keeping
any reaction whose flux clears a cutoff near the solver tolerance. A few
reactions end up carrying flux *just under* that cutoff while doing
load-bearing balancing work, and dropping them leaves a mass-balance residual
of the same order — above the solver's feasibility tolerance. The kept flux
vector is then not actually feasible in the extracted model, and the chains
that relied on those sub-cutoff reactions can carry no flux at all — **a core
reaction included**.

On Recon3D this is measured, not hypothetical: see
:doc:`../tutorials/recon3d`. The mixed-integer formulations should not have
the problem, since a binary at zero forces its reaction's flux to exactly zero
and so discards nothing.

The damaging case is a *blocked core reaction*. It is present in the model, so
its presence alone suggests the extraction worked, but it cannot carry flux
and so cannot play the role it was selected for.
``report.blocked_core`` names them.

Reading the report
------------------

.. code-block:: python

   >>> report
   <ExtractionReport 27 reactions, 0 core missing, 0 blocked (0 of them core), valid=True>

   >>> print(report.summary())
   27 reactions in the extracted model
     every core reaction is present
     flux consistent: every reaction can carry flux

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Attribute
     - Meaning
   * - ``is_valid``
     - Both promises kept. This is the one to assert on.
   * - ``is_consistent``
     - Nothing in the model is blocked.
   * - ``missing_core``
     - Core reactions absent from the model. Should always be empty.
   * - ``blocked``
     - Reactions present but unable to carry flux.
   * - ``blocked_core``
     - The core reactions among them — the damaging ones.
   * - ``n_reactions``
     - Size of the extracted model.

In a pipeline
-------------

.. code-block:: python

   report = check_extraction(extracted, core, tol=1e-4)
   assert report.is_valid, report.summary()

The check solves a handful of LPs, so it is cheap next to the extraction
itself. Run it whenever the result matters.

Reference
---------

.. autofunction:: spectra_cobra.check_extraction

.. autoclass:: spectra_cobra.ExtractionReport
   :members:

See also
--------

* :doc:`../tutorials/recon3d` — the problem measured at genome scale
* :doc:`formulations` — which formulations are immune
