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
any reaction whose flux clears a cutoff near the solver tolerance. An L1
objective has many optimal solutions and spreads tiny fluxes across thousands
of reactions, so discarding everything below that cutoff throws away flux that
was balancing a hub metabolite. The kept flux vector is then only
*approximately* mass balanced, and a reaction whose own flux was near ``tol``
can turn out unable to carry it — **a core reaction included**.

On Recon3D this is not hypothetical: see :doc:`../tutorials/recon3d` for
measured numbers. The mixed-integer formulations do not have the problem,
because a binary at one forces its reaction to carry at least ``tol`` and a
binary at zero forces exactly zero.

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
