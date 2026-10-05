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

Growth is a bound, not a tolerance
----------------------------------

``min_growth`` puts a lower bound on the biomass reaction. ``tol`` stays at
the value any extraction would use.

The alternative — making biomass the core reaction and raising ``tol`` to
the growth rate — appears to work and does not. ``tol`` is the flux *every*
core reaction must carry, so setting it to 0.1 to mean "grow at 0.1" also
demands 0.1 through every other core reaction you name, which is rarely
what you meant. Keeping them apart also makes the answer stop moving: on
iJO1366 the same 15 reactions come back at every tolerance from ``1e-4`` to
``1e-7``, where the conflated version silently returned a model that did
not grow at all below ``0.01``.

``min_growth`` defaults to 0.1, a plausible rate rather than a measured one.
Pass your own where you have one.

Why tradeOff is not offered
---------------------------

:func:`~spectra_cobra.gapfill_for_growth` refuses ``problem_type="tradeOff"``
rather than letting it fail at solve time. ``tradeOff`` requires every
included reaction to carry at least ``tol``, and a growth solution needs
trace fluxes far below that — on iJO1366, 24 of the 442 active reactions
carry less than ``1e-4``, the smallest ``1.96e-06``. It is infeasible at
every tolerance, with or without the growth bound.

Use it for :doc:`../tutorials/gapfilling_tasks`, where the requirement is a
set of reactions rather than biomass, and it works.

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

:doc:`../tutorials/multicellular` works this through, on a toy and then at
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
* :doc:`../tutorials/multicellular` — communities and tissues
* :doc:`community` — building the community first
* :doc:`tasks` — what a task is
