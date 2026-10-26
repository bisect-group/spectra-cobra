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

The nine reconstructions are AGORA 1.03, which is public. This downloads
them, so the page runs from a clean directory with nothing to obtain by
hand:

.. code-block:: python

   from pathlib import Path
   from urllib.request import urlretrieve

   from cobra.io import load_matlab_model
   from spectra_cobra import build_community_model

   AGORA = ("https://raw.githubusercontent.com/VirtualMetabolicHuman/AGORA"
            "/master/CurrentVersion/AGORA_1_03/AGORA_1_03_mat")

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
           urlretrieve(f"{AGORA}/{name}.mat", target)

   models = [load_matlab_model(str(models_dir / f"{name}.mat"))
             for name in strains]
   community = build_community_model(models, organisms=tags)

.. note::

   Building this prints nine warnings, one per organism, saying the anchor
   cannot carry flux and the unit is inert. That is correct and expected:
   AGORA reconstructions arrive with their exchanges closed, so nothing can
   grow until a diet opens them. The next section does that.

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

.. code-block:: text

   2 of 9 organisms; growth 3.694 (needed 2.915), production 8.657 (needed 7.02)
   ['Ecol', 'Erec'] ['Badl', 'Bthe', 'Efae', 'Fpra', 'Kpne', 'Lcas', 'Sthe']

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

   for forced in (["Fpra"], ["Erec"]):
       result = minimal_microbiome(community, products=["EX_but_u"],
                                   required=forced, growth_fraction=0.8,
                                   product_fraction=0.8)
       print(f"required={forced}: {sorted(result.present)}")
       print("   ", result.summary())

.. code-block:: text

   required=['Fpra']: ['Ecol', 'Fpra']
        2 of 9 organisms; growth 3.694 (needed 2.915), production 14.46 (needed 7.02)
   required=['Erec']: ['Ecol', 'Erec']
        2 of 9 organisms; growth 3.694 (needed 2.915), production 8.657 (needed 7.02)

Both are two-member communities meeting every constraint, so both are
optimal: the objective counts organisms, and both cost two. Note that the
published pair makes considerably *more* butyrate — 14.5 against 8.7 — and
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
* Magnúsdóttir, S., Heinken, A., Kutt, L., et al. (2017). Generation of
  genome-scale metabolic reconstructions for 773 members of the human gut
  microbiota. *Nature Biotechnology*, 35(1), 81–89.
  https://doi.org/10.1038/nbt.3703

Next
----

* :doc:`../functions/community` — the API
* :doc:`gapfilling_media` — the other way of saying "the model must work"
