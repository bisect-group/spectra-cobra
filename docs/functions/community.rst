Communities and minimal microbiomes
===================================

:func:`~spectra_cobra.build_community_model` joins organism models into a
community that trades through a shared pool.
:func:`~spectra_cobra.minimal_microbiome` reduces such a community to the
smallest subset that still does the job.

Building a community
--------------------

.. code-block:: python

   from spectra_cobra import build_community_model

   community = build_community_model([model_a, model_b], organisms=["A", "B"])
   community.model            # the joined cobra model
   community.organisms        # ("A", "B")
   community.biomass_reactions  # {"A": "biomass__A", ...}

Each organism keeps its own external compartment; a transport links it to
the shared pool, and one community exchange connects the pool to the
environment. A metabolite therefore moves organism → pool → organism, which
is what makes cross-feeding possible.

Two details carry the weight.

**The organisms' own exchanges are replaced, not kept.** Leaving them in
would let every organism draw on the environment directly, and no amount of
community structure would constrain anything.

**Every reaction is coupled to its organism's biomass**, so an organism that
is not growing carries no flux at all. Without it a dead producer goes on
feeding its neighbours for free, which the test suite checks by building a
community both ways.

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
organisms. The biomass coupling supplies the rest: an organism whose binary
is off cannot grow, and so does nothing.

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

.. autoclass:: spectra_cobra.CommunityModel
   :members:

.. autofunction:: spectra_cobra.minimal_microbiome

.. autoclass:: spectra_cobra.MinimalMicrobiome
   :members:

See also
--------

* :doc:`spectra_me` — the extraction this is built on
* :doc:`formulations` — ``indicator_reactions`` and what it changes
