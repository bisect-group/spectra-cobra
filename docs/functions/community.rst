Communities, tissues and minimal microbiomes
============================================

:func:`~spectra_cobra.build_community_model` joins models into one network
whose parts trade through a pool. The parts may be organisms of a community
or tissues of a body; nothing in the construction distinguishes them.
:func:`~spectra_cobra.build_multi_tissue_model` is the same routine with one
reference network replicated for you.
:func:`~spectra_cobra.minimal_microbiome` reduces a community to the
smallest subset that still does the job.

:doc:`../tutorials/multicellular` builds all of this up from a toy and
explains pools, transports, exchanges and coupling. This page is the
reference.

At a glance
-----------

.. code-block:: python

   from spectra_cobra import build_community_model

   community = build_community_model([model_a, model_b], organisms=["A", "B"])
   community.model               # the joined cobra model
   community.biomass_reactions   # {"A": "biomass__A", ...}
   community.decompose()         # back to one model per unit
   print(community.coupling_summary())

The inputs are read, never modified and never copied.

Choosing the arguments
----------------------

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Argument
     - Reach for it when
   * - ``mode``
     - ``"shared"`` for a large community; ``"pooled"`` when units must
       meet in different places or you have per-unit rates to impose
   * - ``pools``
     - the units do not all meet in one place
   * - ``environment``
     - required as soon as there is more than one pool
   * - ``databases``
     - the units are drafts to be gap-filled against a universal model
   * - ``pool_medium``
     - what the environment supplies
   * - ``link_bounds``
     - you have measured uptake or secretion rates per unit
   * - ``biomass_reactions``
     - the anchor is not the objective and is not named ``biomass*``
   * - ``couple``
     - leave on for organisms; think first for tissues

Finding a minimal microbiome
----------------------------

.. code-block:: python

   from spectra_cobra import minimal_microbiome

   result = minimal_microbiome(community, products=["EX_but_u"])
   result.present      # the organisms kept
   result.membership   # every organism, kept or not

This is the extraction machinery pointed at organisms instead of reactions.
Each organism's biomass reaction gets a binary through
``indicator_reactions`` and a weight of 1, everything else is weighted 0,
and minimising the weighted count becomes minimising the number of
organisms. The coupling supplies the rest: an organism whose binary is off
cannot grow, and so does nothing.

The growth and production requirements are constraints rather than part of
the objective, and they are *measured before they are imposed*: the full
community is solved for its growth and for how much it can produce, and the
minimal community is held to a fraction of each. Asking for a fraction of
something unmeasured would make the answer depend on the units of the
biomass reaction.

.. note::

   The production a result reports is the most the minimal community can
   make while still growing, not whatever it happened to make at its growth
   optimum. The two differ: at the growth optimum the solver has no reason
   to produce anything, so reading production there reports zero for a
   community that in fact meets its requirement.

Reference
---------

.. autofunction:: spectra_cobra.build_community_model

.. autofunction:: spectra_cobra.build_multi_tissue_model

.. autoclass:: spectra_cobra.CommunityModel
   :members:

.. autofunction:: spectra_cobra.minimal_microbiome

.. autoclass:: spectra_cobra.MinimalMicrobiome
   :members:

See also
--------

* :doc:`../tutorials/multicellular` — the worked explanation
* :doc:`gapfilling` — filling the gaps in every unit at once
* :doc:`spectra_me` — the extraction this is built on
* :doc:`formulations` — ``indicator_reactions`` and what it changes
