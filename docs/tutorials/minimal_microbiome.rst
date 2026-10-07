Minimal microbiomes
===================

A gut community contains many organisms that contribute nothing to the
function you care about. A *minimal microbiome* is the smallest subset that
still performs it — still grows, and still makes the metabolite you are
interested in — at a stated fraction of what the whole community managed.

This is the extraction machinery pointed at organisms rather than reactions,
and the mapping is exact. Give each organism's biomass reaction a binary and
a weight of 1, weight everything else 0, and minimising the weighted count
*is* minimising the number of organisms present. It runs through
:func:`~spectra_cobra.spectra_me` with ``problem_type="minNetMILP"`` and
``indicator_reactions`` set to the biomass reactions — the same routine as
every other reconstruction in the package, with a different objective
pointed at a different thing.

The membership comes back through ``return_solutions``, not off the reduced
model. The distinction matters: an organism's biomass reaction can survive
in the extracted model on a trace of flux while its binary is off, so
reading presence from the model would count an absent organism as present.
``MilpSolution.selected`` reports the binaries themselves.

Building a community
--------------------

.. code-block:: python

   from cobra.io import load_matlab_model
   from spectra_cobra import build_community_model

   models = [load_matlab_model(f"agora/{name}.mat") for name in strains]
   community = build_community_model(models, organisms=tags)

Each organism keeps its own external compartment; a transport links it to a
shared pool, and one community exchange connects the pool to the
environment. So a metabolite moves organism → pool → organism, and
cross-feeding becomes possible.

.. warning::

   Check that the pool is actually shared. Collections differ on whether
   the external compartment is written ``glc_D_e`` or ``glc_D[e]``, and if
   the builder fails to recognise one of them each organism ends up
   trading through its own private pool named after its own spelling.
   Nothing errors; the community simply never cross-feeds. Both spellings
   are handled, but it is worth confirming the count of shared exchanges
   looks like a community rather than a sum of isolates.

Applying a diet
---------------

A diet is a set of uptake bounds on the *shared* exchanges:

.. code-block:: python

   for exchange in community.community_exchanges:
       community.model.reactions.get_by_id(exchange).lower_bound = 0.0
   for metabolite, uptake in diet.items():
       rxn_id = f"EX_{metabolite}_u"
       if rxn_id in community.model.reactions:
           community.model.reactions.get_by_id(rxn_id).lower_bound = -uptake

Close everything first. A diet is defined as much by what it withholds as by
what it supplies, and the organisms' default bounds are usually far more
permissive than you intend.

Reducing the community
----------------------

.. code-block:: python

   from spectra_cobra import minimal_microbiome

   result = minimal_microbiome(
       community,
       products=["EX_but_u"],        # the function to preserve
       growth_fraction=0.8,          # keep 80% of the community's growth
       product_fraction=0.8,         # and 80% of its production
   )
   print(result.summary())
   print(sorted(result.present), sorted(result.absent))

The requirements are **measured before they are imposed**: the full
community is solved for its growth, then for how much it can produce while
each organism keeps ``growth_optimum_fraction`` of its own growth, and the
minimal community is held to a fraction of each. Asking for a fraction of
something unmeasured would make the answer depend on the units of the
biomass reaction.

A published comparison
----------------------

The nine-member gut community of Raghu et al. — *Bacteroides
thetaiotaomicron*, *Eubacterium rectale*, *Faecalibacterium prausnitzii*,
*Enterococcus faecalis*, *Lactobacillus casei*, *Streptococcus
thermophilus*, *Bifidobacterium adolescentis*, *Escherichia coli* and
*Klebsiella pneumoniae* — on a Western diet, preserving butyrate.

The community built here has 12162 reactions over 10168 metabolites with
341 shared exchanges. Its behaviour reproduces the published figures:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * -
     - Published
     - Here
   * - community growth (h\ :sup:`-1`)
     - 3.643
     - 3.6433
   * - community butyrate (mmol/gDW-h)
     - 8.775
     - 8.775
   * - 80% growth requirement
     - 2.915
     - 2.915
   * - 80% butyrate requirement
     - 7.020
     - 7.020
   * - minimal microbiome size
     - 2
     - 2

The membership differs: the published pair is *E. coli* with *F.
prausnitzii*, and the pair found here is *E. coli* with *E. rectale*. That
is a tie rather than a disagreement, and it can be checked directly by
forcing each in:

.. code-block:: python

   for forced in (["Fp"], ["Er"]):
       result = minimal_microbiome(community, products=["EX_but_u"],
                                   required=forced, growth_fraction=0.8,
                                   product_fraction=0.8)

.. code-block:: text

   required=['Fp']: ['Ec', 'Fp']  butyrate 13.669 (needed 7.020), growth 3.570
   required=['Er']: ['Ec', 'Er']  butyrate 8.571 (needed 7.020), growth 3.584

Both are two-member communities meeting every constraint, so both are
optimal: the objective counts organisms, and both cost two. Note that the
published pair makes considerably *more* butyrate — 13.7 against 8.6 — and
that this does not break the tie, because production is a constraint here
and not the thing being minimised. If you want the most productive of the
smallest communities, that is a second optimisation, not this one.

.. note::

   Ties are the normal case, not an edge case. A gut community is
   functionally redundant by nature, which is the whole reason a minimal
   microbiome is interesting, and the same redundancy means several
   subsets of equal size usually qualify. Treat a single run as *one*
   minimal microbiome rather than *the* minimal microbiome.

Keeping an organism whatever the objective says
-----------------------------------------------

``required`` forces organisms in. Use it to test a hypothesis — is this
species replaceable? — or to hold a strain you know must be present:

.. code-block:: python

   result = minimal_microbiome(community, products=["EX_but_u"],
                               required=["Fp"])

Next
----

* :doc:`../functions/community` — the API
* :doc:`gapfilling_media` — the other way of saying "the model must work"
