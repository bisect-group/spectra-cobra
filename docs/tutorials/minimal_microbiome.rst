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

The nine reconstructions are **AGORA2 version 2.0**, built 20 June 2022 by
the Thiele group and obtained from `VMH <https://vmh.life/>`_. The copies
used here are the ones committed alongside the SPECTRA manuscript, which
makes them fetchable in a few lines so the page runs from a clean directory
with nothing to obtain by hand:

.. code-block:: python

   from pathlib import Path
   from urllib.request import urlretrieve

   from cobra.io import load_matlab_model
   from spectra_cobra import build_community_model

   AGORA2 = ("https://raw.githubusercontent.com/bisect-group/SPECTRA"
             "/paper_codes/MinimalMicrobiome/models")

   strains = [
       "Bacteroides_thetaiotaomicron_VPI_5482",
       "Bifidobacterium_adolescentis_ATCC_15703",
       "Enterococcus_faecalis_V583",
       "Escherichia_coli_SE11",
       "Eubacterium_rectale_ATCC_33656",
       "Faecalibacterium_prausnitzii_A2_165",
       "Klebsiella_pneumoniae_pneumoniae_MGH78578",
       "Lactobacillus_casei_ATCC_334",
       "Streptococcus_thermophilus_LMG_18311",
   ]
   tags = ["Bthe", "Badl", "Efae", "Ecol", "Erec",
           "Fpra", "Kpne", "Lcas", "Sthe"]

   models_dir = Path("agora_models")
   models_dir.mkdir(exist_ok=True)
   for name in strains:
       target = models_dir / f"{name}.mat"
       if not target.exists():
           urlretrieve(f"{AGORA2}/{name}.mat", target)

   models = [load_matlab_model(str(models_dir / f"{name}.mat"))
             for name in strains]
   community = build_community_model(models, organisms=tags)

.. note::

   Building this prints nine warnings, one per organism, saying the anchor
   cannot carry flux and the unit is inert. That is correct and expected:
   AGORA reconstructions arrive with their exchanges closed, so nothing can
   grow until a diet opens them. The next section does that.

   The result has 16,199 reactions over 13,730 metabolites, with 708 shared
   exchanges.

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

A diet is a set of uptake bounds on the *shared* exchanges. The AGORA diets
are tab-separated ``EX_<metabolite>(e)``/bound pairs. The one used here is
the sample diet published alongside the original implementation at
`RamanLab/minMicrobiome <https://github.com/RamanLab/minMicrobiome>`_, which
is also public:

.. code-block:: python

   DIET_URL = ("https://raw.githubusercontent.com/RamanLab/minMicrobiome"
               "/main/tutorial/sample_diet.txt")

   if not Path("sample_diet.txt").exists():
       urlretrieve(DIET_URL, "sample_diet.txt")

   diet = {}
   with open("sample_diet.txt") as handle:
       for line in handle:
           rxn_id, bound = line.split()
           metabolite = rxn_id[len("EX_"):].removesuffix("(e)")
           diet[metabolite] = abs(float(bound))

   print(len(diet), "metabolites in the diet")

.. code-block:: text

   175 metabolites in the diet

.. code-block:: python

   for exchange in community.community_exchanges:
       community.model.reactions.get_by_id(exchange).lower_bound = 0.0
   for metabolite, uptake in diet.items():
       rxn_id = f"EX_{metabolite}_u"
       if rxn_id in community.model.reactions:
           community.model.reactions.get_by_id(rxn_id).lower_bound = -uptake

Close everything first. A diet is defined as much by what it withholds as by
what it supplies, and the organisms' default bounds are usually far more
permissive than you intend. Of the 175 metabolites, 146 match a shared
exchange; the rest name things none of these nine organisms exchange.

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

.. code-block:: text

   3 of 9 organisms; growth 3.289 (needed 2.631), production 16.24 (needed 9.781)
   ['Ecol', 'Efae', 'Fpra'] ['Badl', 'Bthe', 'Erec', 'Kpne', 'Lcas', 'Sthe']

A warning about reactions carrying less flux than the solver can resolve is
expected here and does not affect the answer. It concerns which *reactions*
the solution puts in the model; the membership you asked for is read from
the organism binaries, which is why ``present`` and ``absent`` are reported
from ``MilpSolution.selected`` rather than from the reactions that survived.

The requirements are **measured before they are imposed**: the full
community is solved for its growth, then for how much it can produce while
each organism keeps ``growth_optimum_fraction`` of its own growth, and the
minimal community is held to a fraction of each. Asking for a fraction of
something unmeasured would make the answer depend on the units of the
biomass reaction.

Against the published result
----------------------------

The nine-member gut community is that of Raghu et al. — *Bacteroides
thetaiotaomicron*, *Eubacterium rectale*, *Faecalibacterium prausnitzii*,
*Enterococcus faecalis*, *Lactobacillus casei*, *Streptococcus
thermophilus*, *Bifidobacterium adolescentis*, *Escherichia coli* and
*Klebsiella pneumoniae* — on a Western diet, preserving butyrate.

The published work used **AGORA 1.03**; this page uses **AGORA2**. Those are
different reconstructions of the same nine strains, so the numbers differ
throughout, and the comparison is only meaningful if the reconstruction is
held fixed:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * -
     - Published
     - AGORA 1.03
     - AGORA2
   * - reactions in the community
     - —
     - 12,162
     - 16,199
   * - shared exchanges
     - —
     - 341
     - 708
   * - community growth (h\ :sup:`-1`)
     - 3.643
     - 3.6434
     - 3.2886
   * - community butyrate (mmol/gDW-h)
     - 8.775
     - 8.7755
     - 12.2264
   * - minimal microbiome size
     - 2
     - **2**
     - **3**

Run on AGORA 1.03 — the inputs the paper used — SPECTRA reproduces the
published figures to four decimals and finds a two-member community. On
AGORA2 the same code finds three. That is the reconstruction talking, not
the algorithm: AGORA2 gives these strains different exchange capabilities,
and the butyrate route that two organisms covered in 1.03 takes three here.

Expect this whenever a reconstruction set is upgraded. If you are comparing
against a published minimal microbiome, match the model version first.

Ties are the normal case
~~~~~~~~~~~~~~~~~~~~~~~~

Several subsets of equal size usually qualify, and which one you get is not
meaningful on its own. Forcing a member in shows the alternatives:

.. code-block:: python

   for forced in (["Fpra"], ["Erec"], ["Ecol"]):
       result = minimal_microbiome(community, products=["EX_but_u"],
                                   required=forced, growth_fraction=0.8,
                                   product_fraction=0.8)
       print(f"required={forced}: {sorted(result.present)}")
       print("   ", result.summary())

.. code-block:: text

   required=['Fpra']: ['Ecol', 'Erec', 'Fpra']
        3 of 9 organisms; growth 3.289 (needed 2.631), production 11.39 (needed 9.072)
   required=['Erec']: ['Ecol', 'Erec', 'Fpra']
        3 of 9 organisms; growth 3.289 (needed 2.631), production 11.39 (needed 9.072)
   required=['Ecol']: ['Ecol', 'Efae', 'Fpra']
        3 of 9 organisms; growth 3.289 (needed 2.631), production 16.24 (needed 9.072)

Two distinct three-member communities, both optimal: the objective counts
organisms and both cost three. They differ considerably in what they make —
16.24 against 11.39 mmol/gDW-h of butyrate — and that does **not** break the
tie, because production is a constraint here and not the thing being
minimised. If you want the most productive of the smallest communities, that
is a second optimisation, not this one.

.. note::

   A gut community is functionally redundant by nature, which is the whole
   reason a minimal microbiome is interesting, and the same redundancy means
   ties are routine. Treat a single run as *one* minimal microbiome rather
   than *the* minimal microbiome.

.. note::

   ``required`` also changes what the requirement *is*. The floors are
   measured on the full community with the forced organisms already pinned,
   so the butyrate requirement above is 9.072 rather than the 9.781 of the
   unforced run. The runs are therefore comparable in size but not in their
   production figures.

Keeping an organism whatever the objective says
-----------------------------------------------

``required`` forces organisms in. Use it to test a hypothesis — is this
species replaceable? — or to hold a strain you know must be present:

.. code-block:: python

   result = minimal_microbiome(community, products=["EX_but_u"],
                               required=["Fpra"])

The organisms are named by the ``tags`` given to
:func:`~spectra_cobra.build_community_model`, so ``required`` takes those
tags and not the strain names. An unknown one is an error rather than a
silent miss:

.. code-block:: text

   SpectraError: Not organisms of this community: ['Fp'].

References
----------

* Raghu, A. K., Palanikumar, I., and Raman, K. (2024). Designing
  function-specific minimal microbiomes from large microbial communities.
  *npj Systems Biology and Applications*, 10, 46.
  https://doi.org/10.1038/s41540-024-00373-1
* Heinken, A., Hertel, J., Acharya, G., et al. (2023). Genome-scale
  metabolic reconstruction of 7,302 human microorganisms for personalized
  medicine. *Nature Biotechnology*, 41(9), 1320–1331.
  https://doi.org/10.1038/s41587-022-01628-0
* Magnúsdóttir, S., Heinken, A., Kutt, L., et al. (2017). Generation of
  genome-scale metabolic reconstructions for 773 members of the human gut
  microbiota. *Nature Biotechnology*, 35(1), 81–89.
  https://doi.org/10.1038/nbt.3703

Next
----

* :doc:`../functions/community` — the API
* :doc:`gapfilling_media` — the other way of saying "the model must work"
