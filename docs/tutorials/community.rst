Microbial community models
==========================

A community is several prokaryotes sharing an environment.

What a community model does is put each member's network side by side in
one model, tag each one so its reactions cannot be confused with another's,
and give them somewhere to pass metabolites to each other. One organism's
secretion arrives there and becomes available to the rest, which is what
makes cross-feeding possible and what makes the community more than the sum
of its members.

This page builds one from nothing on a toy, then gap-fills a 104-organism
community at genome scale.

Organisms, pools and exchanges
------------------------------

Three words carry the whole design.

An **organism** is one member's whole network. Everything belonging to it
gets its identifier as a suffix, so ``PFK`` in organism ``A`` becomes
``PFK__A``. Organisms never share a reaction, even where two of them carry
the same one: *E. coli*'s ``PFK`` and *Bacteroides*' ``PFK`` are separate
reactions with separate fluxes, because the two organisms run them
independently.

A **pool** is a compartment that several organisms can reach. It is the only
way one organism's output becomes another's input. A community usually has a
single pool, the medium everyone sits in, but it may have more, where some
members are in contact with each other and not with the rest, such as a
biofilm layer that only its neighbours reach.

An **exchange** connects a pool to the environment outside the whole system.
Only the pools you nominate get them. This is deliberate and it matters: a
pool that exists purely to link two organisms must stay closed, or those two
could dump metabolites straight out of the system and bypass everything
else.

Every diagram on this page uses the same five shapes:

.. graphviz::
   :align: center

   digraph legend {
     rankdir=LR; bgcolor="transparent"; nodesep=0.18; ranksep=0.30;
     node [fontname="Helvetica", fontsize=9.5];
     edge [style=invis];

     met  [label="metabolite\nin an organism", shape=ellipse, style=filled,
           fillcolor="#e8f0fa", color="#3c6997", fontcolor="#1b3a57"];
     pool [label="pooled\nmetabolite", shape=ellipse, style=filled,
           fillcolor="#fdf0e3", color="#c07d2a", fontcolor="#6b4411",
           penwidth=1.6];
     rxn  [label="reaction", shape=box, style="filled,rounded",
           fillcolor="#ffffff", color="#3c6997", fontcolor="#1b3a57",
           height=0.34];
     gone [label="reaction a\ndraft lacks", shape=box,
           style="filled,rounded,dashed", fillcolor="#fdf2ef", color="#b3411f",
           fontcolor="#b3411f", height=0.34];
     organism [label="an organism", shape=box, style="rounded,filled",
           fillcolor="#f6f9fd", color="#c7d6e8", fontcolor="#3c6997",
           height=0.34];

     met -> pool -> rxn -> gone -> organism;
   }

Blue is private to an organism, orange is shared, and a dashed red outline marks
something a draft does not have.

Note what has to *go*: each organism's own exchange reactions are removed. A
draft model has ``EX_glc_e`` so it can eat; leave that in place and every
organism feeds itself directly from the environment, the pool is decorative, and
nothing is a community at all.

How the pool is built
---------------------

``mode="shared"`` is what a community wants: the organisms share the external
compartment itself. One row of the stoichiometric matrix per exchangeable
metabolite, for the whole community, and no transport reactions at all.

.. graphviz::
   :caption: ``mode="shared"``. There is one ``glc[e]`` and both organisms
             reach it directly. Smaller, but with no transport between organism
             and pool there is nothing for ``link_bounds`` to attach a
             per-organism uptake rate to.
   :align: center

   digraph shared {
     rankdir=LR; bgcolor="transparent"; compound=true; nodesep=0.28; ranksep=0.5;
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#6b7280", arrowsize=0.7];

     subgraph cluster_a {
       label="organism A"; labeljust=l; fontname="Helvetica"; fontsize=9;
       fontcolor="#3c6997"; color="#c7d6e8"; style="rounded,filled";
       fillcolor="#f6f9fd"; margin=10;
       upA [label="UP", shape=box, style="filled,rounded", width=0.7,
            height=0.26, fillcolor="#ffffff", color="#3c6997",
            fontcolor="#1b3a57"];
     }
     subgraph cluster_b {
       label="organism B"; labeljust=l; fontname="Helvetica"; fontsize=9;
       fontcolor="#3c6997"; color="#c7d6e8"; style="rounded,filled";
       fillcolor="#f6f9fd"; margin=10;
       upB [label="UP", shape=box, style="filled,rounded", width=0.7,
            height=0.26, fillcolor="#ffffff", color="#3c6997",
            fontcolor="#1b3a57"];
     }

     pool [label="glc[e]", shape=ellipse, style=filled, fillcolor="#fdf0e3",
           color="#c07d2a", fontcolor="#6b4411", penwidth=1.6];
     ex   [label="EX_glc[e]", shape=box, style="filled,rounded", width=0.9,
           height=0.26, fillcolor="#ffffff", color="#c07d2a",
           fontcolor="#6b4411"];
     env  [label="environment", shape=plaintext, fontcolor="#9ca3af", fontsize=9];

     upA -> pool [dir=both, color="#c07d2a"];
     upB -> pool [dir=both, color="#c07d2a"];
     pool -> ex [dir=both, color="#c07d2a"];
     ex -> env [dir=both, style=dashed, color="#c9ced6"];
   }

.. note::

   The other mode, ``mode="pooled"``, gives each organism a private external
   compartment and a transport reaction into the pool. That buys two things a
   community rarely needs and always pays for: **several pools**, so a member
   can meet one group of organisms in one place and another group somewhere
   else, and **per-organism bounds** on what crosses into each pool, through
   ``link_bounds``.

   It costs one transport per organism per exchangeable metabolite. A hundred
   organisms across a thousand metabolites is a hundred thousand extra
   reactions in ``"pooled"`` and none in ``"shared"`` — the difference
   between a model you can build and one you cannot.

   :doc:`multi_tissue` is where that structure earns its keep, and shows it
   drawn out: tissues genuinely do not all meet in the same place, and the
   uptake rates you measure are per tissue.

Which metabolites are exchangeable
----------------------------------

Reconstructions do not agree on how to write an extracellular metabolite.
BiGG says ``glc__D_e``, VMH and AGORA say ``glc_D[e]``, Human-GEM says
``MAM01965e``. Two questions follow, and they have different answers.

**Which compartment is the outside one?** Worked out per organism, by cobrapy's
own heuristic: a compartment with a recognised name, else the one carrying
the most boundary reactions. Per organism, not once for the community, because
the organisms need not agree — one model may call it ``e`` and another
``extracellular``. Override it when you need to:

.. code-block:: text

   build_community_model(..., external_compartment="ext")          # all organisms
   build_community_model(..., external_compartment={"A": "e", "B": "ext"})

This step keys off ``met.compartment``, not off the identifier. A model
whose metabolites carry no compartment at all — some ``.mat`` dumps — will
have nothing recognised as external, so set them before you build.

**Which metabolites are the same compound?** Two organisms meet in the pool
only if their metabolites resolve to the same name. By default the
compartment is stripped, which resolves ``glc_D_e`` and ``glc_D[e]`` to the
same ``glc_D``, so the two spellings share one pool row in either mode.

What it cannot do is reconcile identifiers from *different namespaces*.
``glc__D_e`` and ``glc_D[e]`` are the same molecule under two schemes, and
no amount of suffix-stripping will tell you so. Left alone that failure is
silent and total: every organism trades through a pool of its own and nothing
cross-feeds, with no error anywhere. So the builder checks, and says so:

.. code-block:: text

   WARNING pool 'u' is shared by 2 organisms but not one metabolite in it is
   reached by more than one of them, so nothing can cross-feed. The usual
   cause is identifiers from different namespaces; pass metabolite_key to
   say how they correspond. Examples: {'A': ['glc__D_e'], 'B': ['glc_D[e]']}

``metabolite_key`` is the way out. Given a metabolite, return the name of
the thing it *is*:

.. code-block:: python

   import re
   from cobra.io import load_matlab_model, load_model

   # agora.mat stands in for any model of yours in the other namespace;
   # substitute your own file here.
   bigg_model = load_model("textbook")              # spells it ``glc__D_e``
   agora_model = load_matlab_model("agora.mat")     # spells it ``glc_D[e]``

   build_community_model(
       [bigg_model, agora_model],
       organisms=["A", "B"],
       metabolite_key=lambda met: re.sub(r"_+D(_e|\[e\])$", "_D", met.id),
   )

A mapping table from one namespace to the other does the same job and is
what you want for anything real.

.. tip::

   ``pool_medium`` and ``link_bounds`` are keyed by metabolite identifier
   as *your* models spell it, and are resolved by looking each one up
   among the metabolites that exist rather than by stripping suffixes. A
   name matching nothing is reported rather than silently ignored, which
   is how you find out a medium file and a model disagree.

A toy community, end to end
---------------------------

A four-step chain, ``glc → a → b → c``, with an anchor that consumes ``c``.
Two drafts each hold part of it:

.. graphviz::
   :caption: The same chain seen twice, once per draft. Dashed red is what
             that draft lacks: ``A`` gets as far as ``b`` but cannot make
             ``c``, and ``B`` can finish but cannot start. The pool ``b_u``
             is how ``A`` hands ``b`` to ``B``.
   :align: center

   digraph chain {
     rankdir=LR; bgcolor="transparent"; nodesep=0.22; ranksep=0.34;
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#6b7280", arrowsize=0.6];

     node [shape=ellipse, style=filled, fillcolor="#e8f0fa", color="#3c6997",
           fontcolor="#1b3a57", width=0.34, height=0.30];
     node [shape=box, style="filled,rounded", width=0.56, height=0.26,
           fillcolor="#ffffff", color="#3c6997", fontcolor="#1b3a57"];

     subgraph cluster_a {
       label="draft A — starts, cannot finish"; labeljust=l;
       fontname="Helvetica"; fontsize=9; fontcolor="#3c6997";
       color="#c7d6e8"; style="rounded,filled"; fillcolor="#f6f9fd"; margin=10;
       gA [label="glc", shape=ellipse, fillcolor="#e8f0fa"];
       aA [label="a", shape=ellipse, fillcolor="#e8f0fa"];
       bA [label="b", shape=ellipse, fillcolor="#e8f0fa"];
       cA [label="c", shape=ellipse, fillcolor="#e8f0fa"];
       UPa [label="UP"]; R1a [label="R1"]; BIOa [label="BIO"];
       R2a [label="R2", style="filled,rounded,dashed", color="#b3411f",
            fontcolor="#b3411f", fillcolor="#fdf2ef"];
       gA -> UPa -> aA -> R1a -> bA;
       bA -> R2a [style=dashed, color="#d8b2a6"];
       R2a -> cA [style=dashed, color="#d8b2a6"];
       cA -> BIOa;
     }
     subgraph cluster_b {
       label="draft B — finishes, cannot start"; labeljust=l;
       fontname="Helvetica"; fontsize=9; fontcolor="#3c6997";
       color="#c7d6e8"; style="rounded,filled"; fillcolor="#f6f9fd"; margin=10;
       gB [label="glc", shape=ellipse, fillcolor="#e8f0fa"];
       aB [label="a", shape=ellipse, fillcolor="#e8f0fa"];
       bB [label="b", shape=ellipse, fillcolor="#e8f0fa"];
       cB [label="c", shape=ellipse, fillcolor="#e8f0fa"];
       R2b [label="R2"]; BIOb [label="BIO"];
       UPb [label="UP", style="filled,rounded,dashed", color="#b3411f",
            fontcolor="#b3411f", fillcolor="#fdf2ef"];
       R1b [label="R1", style="filled,rounded,dashed", color="#b3411f",
            fontcolor="#b3411f", fillcolor="#fdf2ef"];
       gB -> UPb [style=dashed, color="#d8b2a6"];
       UPb -> aB [style=dashed, color="#d8b2a6"];
       aB -> R1b [style=dashed, color="#d8b2a6"];
       R1b -> bB [style=dashed, color="#d8b2a6"];
       bB -> R2b -> cB -> BIOb;
     }

     pool [label="b_u", shape=ellipse, style=filled, fillcolor="#fdf0e3",
           color="#c07d2a", fontcolor="#6b4411", penwidth=1.6];
     bA -> pool [dir=both, color="#c07d2a", constraint=false];
     pool -> bB [dir=both, color="#c07d2a", constraint=false];
   }

Build the community. Each organism is backed by a *database* — everything it is
allowed to borrow — and its own draft supplies the bounds:

.. code-block:: python

   from cobra import Metabolite, Model, Reaction
   from spectra_cobra import build_community_model

   RECIPES = {
       "UP": {"glc_e": -1.0, "a_c": 1.0},
       "R1": {"a_c": -1.0, "b_c": 1.0},
       "R2": {"b_c": -1.0, "c_c": 1.0},
       "BIO": {"c_c": -1.0},
   }

   def chain_model(model_id, reaction_ids):
       """The glc -> a -> b -> c chain, holding only the named reactions."""
       model = Model(model_id)
       mets = {m: Metabolite(m, compartment=m.rsplit("_", 1)[1])
               for m in ("glc_e", "a_c", "b_c", "c_c")}
       for rxn_id in reaction_ids:
           reaction = Reaction(rxn_id, lower_bound=0.0, upper_bound=1000.0)
           reaction.add_metabolites(
               {mets[m]: coeff for m, coeff in RECIPES[rxn_id].items()})
           model.add_reactions([reaction])
       uptake = Reaction("EX_glc_e", lower_bound=-10.0, upper_bound=1000.0)
       uptake.add_metabolites({mets["glc_e"]: -1.0})
       model.add_reactions([uptake])
       model.objective = "BIO"
       return model

   chain = chain_model("chain", ["UP", "R1", "R2", "BIO"])   # the database
   draft_a = chain_model("A", ["UP", "R1", "BIO"])           # cannot finish
   draft_b = chain_model("B", ["R2", "BIO"])                 # cannot start

   community = build_community_model(
       [draft_a, draft_b],
       organisms=["A", "B"],
       databases={"A": chain, "B": chain},
       pool_medium={"glc_e": (-10.0, 1000.0)},
   )

   community.database_reactions["A"]   # ('R2__A',)
   community.database_reactions["B"]   # ('R1__B', 'UP__B')

``pool_medium`` is what the environment offers, keyed by metabolite as your
own models spell it. Everything unnamed stays shut to uptake.

Now fill the gaps in both at once:

.. code-block:: python

   from spectra_cobra import gapfill_community

   result = gapfill_community(community, problem_type="minNetMILP")

   result.added        # {'A': ('R2',), 'B': ()}
   result.models()     # one untagged cobra model per organism

``B`` borrows nothing. What it was missing arrives from ``A`` through the
pool. ``A`` still needs its one reaction, because it has to feed itself
whatever ``B`` does. Filling the two separately costs three reactions
instead of one — and that gap is the entire argument for doing it this way.

What the gap-filler charges for
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 44 22 34

   * - Reaction
     - Weight
     - Why
   * - already in the organism's draft
     - 0
     - it is not an addition
   * - only in the organism's database
     - 1
     - this is what is counted
   * - transport into a pool
     - 0
     - free in ``"shared"`` mode too
   * - community exchange
     - depends
     - see below

The exchange weight is the interesting one, and it is the one place where
the two formulations must be told apart.

Under ``minNetMILP`` the objective counts reactions, so charging 1 for an
exchange makes the solver reach for a neighbour's secretion before it
reaches for fresh material from the medium — exactly the behaviour a
community model exists to capture.

Under ``minNetLP`` the objective sums weighted *flux*, and a shared
exchange carries the flux of every organism drawing on it. Charging for it
therefore penalises an organism for having company, and the solver responds by
adding internal reactions so that it can take up less. Measured on hCom
organisms against filling each of them on its own:

.. list-table:: Reactions added, ``minNetLP``
   :header-rows: 1
   :widths: 22 26 26 26

   * - Organisms
     - weight 1
     - weight 0 (default)
     - one at a time
   * - 2
     - 61
     - **8**
     - 9
   * - 4
     - 89
     - **17**
     - 18
   * - 6
     - 21
     - **21**
     - 31

At weight 0 the community beats filling the members separately at every
size, and the sequence is orderly. At weight 1 it is erratic, and only
catches up once there are enough organisms to share the fixed cost of the
exchanges between them — the pool has one set of them however many organisms
there are, so the per-organism charge falls as the community grows. It is also
slower while it is wrong: 2,931 seconds against 854 at six organisms.

The mixed-integer formulation does not have the problem, because its
objective counts reactions rather than flux, so an exchange costs the same
whether one organism uses it or twenty. At two organisms it adds 8 against 10
filled separately, with the weight left at 1.

So the default follows the formulation: 1 for the MILP, 0 for the LP. Pass
``exchange_weight`` to override it either way. The published pipeline uses
1 throughout, and at the scale it reports on — a hundred and four
organisms, with the mixed-integer formulation — that is the right number.
It is small communities under the LP where it misleads.

Two more defaults worth knowing. Each organism's anchor reaction is **core**, so
every organism has to work — otherwise the cheapest community is one where half
the members are dead. And ``keep_draft=True`` keeps every reaction an organism
already had: the formulation returns the *smallest network* meeting the
requirement, which left to itself would throw away parts of the draft that
carry no flux under this particular medium. Gap-filling adds; it should not
quietly subtract.

Coupling, and what an organism is gated on
------------------------------------------

An absent organism must be absent in full. Not merely not growing — carrying no
flux at all. Otherwise a dead organism goes on running its metabolism and
feeding its neighbours for free.

That is what coupling does. For every reaction of an organism,

.. math::

   |v_i| \le c \cdot v_{\text{anchor}} + u

so a reaction can carry flux only in proportion to its organism's **anchor**.
The defaults, :math:`c = 1000` and :math:`u = 0.01`, are the published
values for a microbial community whose anchor is biomass.

Which reaction is the anchor is decided in this order, and never guessed
past it:

1. what you passed in ``biomass_reactions``,
2. the model's objective, if exactly one reaction is in it,
3. a single reaction named ``biomass*``,
4. otherwise it refuses, and tells you what the candidates looked like.

Refusing matters. Coupling a few thousand reactions to the wrong anchor
would throttle or free that whole organism without raising anything, and a
draft reconstruction may carry several biomass-like reactions or none that
is named recognisably.

Whatever is chosen is reported:

.. code-block:: python

   print(community.coupling_summary())

.. code-block:: text

   organism    anchor                  how          capacity        cap
   tis1    DM_atp_c___tis1         given             10.0      1e+04
   tis2    DM_atp_c___tis2         given           0.0001       0.11

``capacity`` is the most that anchor can carry; ``cap`` is what the coupling
then permits every reaction in that organism. The second row is the failure
to watch for. An anchor running at 1e-4 with :math:`c = 1000` caps the whole
organism at 0.11, which strangles it silently. The builder warns when it
sees this.

.. note::

   :math:`c` is not a number to be derived; it is a modelling choice about
   how much flux an organism may carry per unit of anchor flux. If your anchor runs at
   a different scale from a bacterial growth rate, raise :math:`c` to match
   it, or pass ``couple=False`` and do without.

Community gap-filling at genome scale
-------------------------------------

The published case study gap-fills the hCom synthetic gut community against
CarveMe's universal reconstructions. Each organism's compartment holds an
entire universal — the Gram-positive, Gram-negative or common one, by Gram
stain — bounded by its own draft wherever the two agree. The organisms share
the extracellular compartment, constrained by the standard amino acid
complete (SAAC) medium.

.. note::

   This section needs the case study's inputs, which are not part of this
   package. They are on the ``paper_codes`` branch of the `SPECTRA
   repository <https://github.com/bisect-group/SPECTRA/tree/paper_codes>`_:

   * ``hCOM_community/CommunityGapFilling/hComDraftModels/`` — 104 CarveMe
     drafts, one per organism
   * ``hCOM_community/CommunityGapFilling/Universal_models/`` — CarveMe's
     universals, general (``carvemeModelConv``), Gram-positive (``GP``) and
     Gram-negative (``GN``)
   * ``hCOM_community/CommunityGapFilling/gram_stain.csv`` — the stain per
     strain, as ``0``, ``+`` or ``-``
   * ``hCOM_community/SAACmedia.xlsx`` — the medium

   Paths below assume you are in ``hCOM_community/CommunityGapFilling``.

The drafts are CarveMe reconstructions and the universals CarveMe's own,
picked per organism by Gram stain — three of them, because a strain whose
stain is unrecorded gets the general universal rather than a guess:

.. code-block:: python

   from pathlib import Path

   import pandas as pd
   from cobra.io import load_matlab_model

   names = sorted(path.stem for path in Path("hComDraftModels").glob("*.mat"))
   drafts = [load_matlab_model(f"hComDraftModels/{name}.mat") for name in names]

   # '0' is "not recorded", and takes the general universal.
   UNIVERSAL_FILE = {"0": "carvemeModelConv_with_biomass.mat",
                     "+": "carvemeModelGPconv_with_biomass.mat",
                     "-": "carvemeModelGNconv_with_biomass.mat"}
   universals = {stain: load_matlab_model(f"Universal_models/{filename}")
                 for stain, filename in UNIVERSAL_FILE.items()}

   stains = pd.read_csv("gram_stain.csv").set_index("Strain")["stain"]
   universal_for = {name: universals[stains[name]] for name in names}

   # The SAAC medium, as uptake bounds on the shared compartment.
   saac = {row.metabolite: (-row.uptake, 1000.0)
           for row in pd.read_excel("../SAACmedia.xlsx").itertuples()}

.. code-block:: python

   community = build_community_model(
       drafts,
       organisms=names,
       mode="shared",
       databases={name: universal_for[name] for name in names},
       pool_medium=saac,
       couple=False,
   )

   result = gapfill_community(community, tol=1e-5, problem_type="minNetLP")

``couple=False`` because the requirement here is carried by the biomass
bound, not by the coupling: the published pipeline floors biomass and ATP
maintenance at 0.1 so that both are kept, and makes biomass the only core
reaction.

.. note::

   One practical snag, and it is not SPECTRA's. The CarveMe universal
   contains a reaction named ``St`` (sulfur diffusion), and ``St`` is the LP
   file format's keyword for the constraint section. optlang copies a Gurobi
   model by writing an LP file and reading it back, so a model holding that
   reaction produces a file Gurobi cannot parse.

   SPECTRA handles it: when a solver cannot serialise itself the model is
   rebuilt from its dictionary form instead, carrying reactions,
   metabolites, genes, bounds and GPRs but not the solver, which is then
   reattached. You will see a warning saying so. Renaming the reaction
   avoids the detour.

References
----------

* Machado, D., Andrejev, S., Tramontano, M., and Patil, K. R. (2018). Fast
  automated reconstruction of genome-scale metabolic models for microbial
  species and communities. *Nucleic Acids Research*, 46(15), 7542–7553.
  https://doi.org/10.1093/nar/gky537
* Magnúsdóttir, S., Heinken, A., Kutt, L., et al. (2017). Generation of
  genome-scale metabolic reconstructions for 773 members of the human gut
  microbiota. *Nature Biotechnology*, 35(1), 81–89.
  https://doi.org/10.1038/nbt.3703
* Cheng, A. G., Ho, P.-Y., Aranda-Díaz, A., et al. (2022). Design,
  construction, and in vivo augmentation of a complex gut microbiome. *Cell*,
  185(19), 3617–3636. https://doi.org/10.1016/j.cell.2022.08.003

Reference
---------

* :doc:`multi_tissue` — the same construction, applied to tissues
* :doc:`../functions/community` — the builder and the minimal microbiome
* :doc:`../functions/gapfilling` — ``gapfill_community`` and its siblings
* :doc:`minimal_microbiome` — reducing a community instead of filling it
