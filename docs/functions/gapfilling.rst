Gap-filling
===========

A draft reconstruction usually cannot do everything the organism can.
Gap-filling borrows the smallest useful set of reactions from a universal
model to close the difference, and what "the difference" means is the choice
you make.

:func:`~spectra_cobra.gapfill_for_growth` requires the model to grow, in
each of a panel of media. :func:`~spectra_cobra.gapfill_for_tasks` requires
it to perform a list of :doc:`tasks`.

Usage
-----

.. code-block:: python

   from spectra_cobra import gapfill_for_growth

   M9 = {"EX_pi_e": -1000, "EX_h2o_e": -1000, "EX_h_e": -1000,
         "EX_nh4_e": -1000, "EX_o2_e": -1000, "EX_co2_e": -1000}
   media = {carbon: {**M9, f"EX_{carbon}_e": -10.0}
            for carbon in ("glc__D", "ac", "succ", "pyr")}

   result = gapfill_for_growth(draft, universal, media)
   print(result.summary())
   print(result.added)

A medium entry may be a number or a ``(lower, upper)`` pair, and the two
say different things. Uptake is negative flux, so a bare ``-10`` is a
lower bound — the model *may* take up to 10, and may take nothing. A pair
pins both ends, and only a negative upper bound can oblige the model to
consume:

.. code-block:: python

   medium = {
       "EX_glc__D_e": -10.0,          # may take up to 10
       "EX_o2_e": (-20.0, -5.0),      # must take between 5 and 20
   }

Either way a medium governs uptake only. Exchanges it does not name have
their lower bound closed; no exchange has its upper bound touched unless
you give it a pair, so secretion stays as the reconstruction has it.

Filling a whole community at once
---------------------------------

:func:`~spectra_cobra.gapfill_community` fills every unit of a community in
one solve. A gap in one unit can then be closed by another's secretion
instead of by a new reaction, which is why the joint answer is the smaller
one.

.. code-block:: python

   from spectra_cobra import build_community_model, gapfill_community

   community = build_community_model(
       drafts, organisms=names, mode="shared",
       databases={name: universal for name in names}, pool_medium=medium,
   )
   result = gapfill_community(community, problem_type="minNetLP")

   result.added      # the reactions each unit took, by unit
   result.models()   # one untagged, gap-filled model per unit

What a unit already has costs nothing, what only its database has costs
one, and a community exchange costs one as well — so the solver reaches for
a neighbour's secretion before it reaches for the medium. Each unit's
anchor is core, which is what stops the cheapest answer being a community
where half the members are dead.

``keep_draft`` is on by default and is worth understanding: the formulation
returns the *smallest network* that meets the requirement, so left alone it
would discard parts of the draft that carry no flux under this medium.
Gap-filling adds; it should not quietly subtract.

:doc:`../tutorials/community` works this through, on a toy and then at
genome scale.

Reference
---------

.. autofunction:: spectra_cobra.gapfill_for_growth

.. autofunction:: spectra_cobra.gapfill_for_tasks

.. autofunction:: spectra_cobra.gapfill_community

.. autoclass:: spectra_cobra.GapfillResult
   :members:

.. autoclass:: spectra_cobra.CommunityGapfillResult
   :members:

See also
--------

* :doc:`../tutorials/gapfilling_media` — growth in a panel of media
* :doc:`../tutorials/gapfilling_tasks` — towards a task list
* :doc:`../tutorials/community` — communities
* :doc:`../tutorials/multi_tissue` — tissues
* :doc:`community` — building the community first
* :doc:`tasks` — what a task is
