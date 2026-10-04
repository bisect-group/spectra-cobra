SPECTRA for cobrapy
===================

**S**\ calable **P**\ latform for **E**\ xtracting **C**\ onstraint-based
**T**\ op-down **R**\ econstructions and **A**\ nalysis.

``spectra-cobra`` reconstructs metabolic networks from multi-omics data at a
range of biological scales: minimal reactomes, context-specific models,
gap-filled reconstructions, minimal microbiomes, microbial community models
and multi-tissue models. What changes between them is the universal model,
the evidence supplied and the objective chosen, not the routine called. It is
built on `cobrapy <https://github.com/opencobra/cobrapy>`_.

At a glance
-----------

.. code-block:: python

   from cobra.io import load_model
   from spectra_cobra import spectra_cc, spectra_me

   model = load_model("textbook")

   # 1. Drop the blocked reactions.
   consistent = spectra_cc(model, tol=1e-3)

   # 2. Extract a model around the reactions you care about.
   extracted = spectra_me(consistent, core_reactions=["PGI", "PFK", "FBA"],
                          tol=1e-3)

What is in the package
----------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Function
     - Does
   * - :doc:`functions/spectra_cc`
     - Removes the blocked reactions of a model
   * - :doc:`functions/spectra_me`
     - Extracts a context-specific model around a core set
   * - :doc:`functions/spectra_ccme`
     - Both at once, for a model that may be inconsistent
   * - :doc:`functions/formulations`
     - The network inference objectives the extraction can use
   * - :doc:`functions/flux_reducer`
     - Finds a sparse flux distribution
   * - :doc:`functions/verification`
     - Checks an extracted model is consistent and complete
   * - :doc:`functions/tasks`
     - States what a network must be able to do, and finds the reactions
       it cannot do it without

Which function do I want?
-------------------------

* To find out which reactions of a universal model can carry flux at all, use
  :doc:`functions/spectra_cc`.
* To build a context-specific model and your universal model is already flux
  consistent, use :doc:`functions/spectra_me`.
* To build one and you are not sure the universal model is consistent, use
  :doc:`functions/spectra_ccme`. It reports the core reactions it had to drop
  instead of failing.
* To fill the gaps in a draft reconstruction, so that it can perform a list
  of metabolic tasks or grow in a panel of media, see
  :doc:`tutorials/gapfilling_tasks` and :doc:`tutorials/gapfilling_media`.

Citation
--------

   S, P. K., Sridhar, S., Alsmadi, N., Mahadevan, R., & Bhatt, N. P. (2026).
   *Generalist method to reconstruct metabolic networks from multi-omics data
   at large-scale.* bioRxiv. https://doi.org/10.64898/2026.04.02.716249

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation

.. toctree::
   :maxdepth: 2
   :caption: Functions

   functions/spectra_cc
   functions/spectra_me
   functions/spectra_ccme
   functions/formulations
   functions/flux_reducer
   functions/verification
   functions/tasks

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/consistency
   tutorials/extraction
   tutorials/formulations
   tutorials/alternative_solutions
   tutorials/topology
   tutorials/gapfilling_tasks
   tutorials/gapfilling_media
   tutorials/recon3d

.. toctree::
   :maxdepth: 1
   :caption: About

   api
