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
any reaction whose flux clears ``inclusion_cutoff``. That cutoff defaults to
``tol * 1e-7``, deliberately far below the solver's own tolerance: raising it
discards reactions that carry almost no flux but are load-bearing for a mass
balance, and the chains that relied on them then collapse, leaving **core
reactions present in the model but unable to carry flux**.

Keeping the cutoff that low fixes the core reactions at the price of also
keeping some reactions that are there only because of numerical noise, so the
model can still contain a few blocked non-core reactions. That trade is
intentional — see :doc:`../tutorials/recon3d` for the measured effect — and
``check_extraction`` is how you see which side of it you landed on.

The mixed-integer formulations are unaffected either way, since they read
their answer off their binaries rather than off the flux.

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
