Multicellular and community models
==================================

A microbial community and a multi-tissue body are the same construction.
Several copies of a metabolic network sit side by side, each one tagged so
its reactions cannot be confused with another's, and they meet somewhere
they can pass metabolites to each other. What differs between a gut
community and a human body is which copies meet where, and what the meeting
place is called.

This page builds one from nothing on a toy, then uses the same routine twice
at scale: once for a multi-tissue model and once to gap-fill a microbial
community.

Units, pools and exchanges
--------------------------

Three words carry the whole design.

A **unit** is one copy of a network: an organism, or a tissue. Everything
belonging to it gets its identifier as a suffix, so ``PFK`` in unit ``A``
becomes ``PFK__A``. Units never share a reaction.

A **pool** is a compartment that several units can reach. It is the only
way one unit's output becomes another's input. A community usually has a
single pool; a body has several, because the liver and the gut meet the
blood but the ovary and the granulosa cells also meet each other.

An **exchange** connects a pool to the environment outside the whole system.
Only the pools you nominate get them. This is deliberate and it matters: a
pool that exists purely to link two tissues must stay closed, or those two
could dump metabolites straight out of the system and bypass everything
else.

Note what has to *go*: each unit's own exchange reactions are removed. A
draft model has ``EX_glc_e`` so it can eat; leave that in place and every
unit feeds itself directly from the environment, the pool is decorative, and
nothing is a community at all.

Two ways to build a pool
------------------------

``mode="pooled"`` — each unit keeps a private external compartment, and a
transport reaction carries each metabolite between it and the pool:

.. graphviz::
   :caption: ``mode="pooled"``. Each unit's own ``glc[e]`` is distinct, and
             ``TR_`` reactions move glucose into the shared ``glc_u``. The
             transports belong to a unit, so they can be bounded per unit.
   :align: center

   digraph pooled {
     rankdir=LR;
     bgcolor="transparent";
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#555555", arrowsize=0.7];
     compound=true;

     node [shape=ellipse, style=filled, fillcolor="#e8f0fa",
           color="#3c6997", fontcolor="#1b3a57"];
     glcA [label="glc[e]__A"];
     glcB [label="glc[e]__B"];
     node [fillcolor="#fdf0e3", color="#c07d2a", fontcolor="#6b4411"];
     pool [label="glc_u"];

     node [shape=box, width=0.32, height=0.24, style=filled,
           fillcolor="#f7f7f7", color="#999999", fontcolor="#333333"];
     trA [label="TR_glc_u__A"];
     trB [label="TR_glc_u__B"];
     ex  [label="EX_glc_u"];

     node [shape=point, width=0.05, color="#bbbbbb"];
     env;

     subgraph cluster_a {
       label="unit A"; color="#cccccc"; fontname="Helvetica"; fontsize=9;
       glcA;
     }
     subgraph cluster_b {
       label="unit B"; color="#cccccc"; fontname="Helvetica"; fontsize=9;
       glcB;
     }

     glcA -> trA [dir=both];
     trA  -> pool [dir=both];
     glcB -> trB [dir=both];
     trB  -> pool [dir=both];
     pool -> ex [dir=both];
     ex -> env [style=dashed, color="#bbbbbb"];
   }

``mode="shared"`` — the units share the external compartment itself. One row
of the stoichiometric matrix per exchangeable metabolite, for the whole
community, and no transports at all:

.. graphviz::
   :caption: ``mode="shared"``. There is one ``glc[e]`` and both units use
             it. Smaller, but there is no per-unit reaction left to bound.
   :align: center

   digraph shared {
     rankdir=LR;
     bgcolor="transparent";
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#555555", arrowsize=0.7];

     node [shape=box, width=0.32, height=0.24, style=filled,
           fillcolor="#f7f7f7", color="#999999", fontcolor="#333333"];
     upA [label="UP__A"]; upB [label="UP__B"]; ex [label="EX_glc[e]"];

     node [shape=ellipse, style=filled, fillcolor="#fdf0e3",
           color="#c07d2a", fontcolor="#6b4411"];
     pool [label="glc[e]"];

     node [shape=point, width=0.05, color="#bbbbbb"];
     env;

     pool -> upA [dir=both];
     pool -> upB [dir=both];
     pool -> ex [dir=both];
     ex -> env [style=dashed, color="#bbbbbb"];
   }

Which to use:

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * -
     - ``"pooled"``
     - ``"shared"``
   * - pools
     - as many as you like
     - exactly one
   * - per-unit bounds on what crosses
     - yes, through ``link_bounds``
     - no reaction to bound
   * - size
     - one transport per unit per metabolite
     - nothing extra
   * - use it for
     - tissues, measured uptake rates
     - large communities

The size difference is not cosmetic. A hundred organisms across a thousand
exchangeable metabolites is a hundred thousand transport reactions in
``"pooled"`` and none in ``"shared"``.

Which metabolites are exchangeable
----------------------------------

Reconstructions do not agree on how to write an extracellular metabolite.
BiGG says ``glc__D_e``, VMH and AGORA say ``glc_D[e]``, Human-GEM says
``MAM01965e``. Two questions follow, and they have different answers.

**Which compartment is the outside one?** Worked out per unit, by cobrapy's
own heuristic: a compartment with a recognised name, else the one carrying
the most boundary reactions. Per unit, not once for the community, because
the units need not agree — one model may call it ``e`` and another
``extracellular``. Override it when you need to:

.. code-block:: python

   build_community_model(..., external_compartment="ext")          # all units
   build_community_model(..., external_compartment={"A": "e", "B": "ext"})

This step keys off ``met.compartment``, not off the identifier. A model
whose metabolites carry no compartment at all — some ``.mat`` dumps — will
have nothing recognised as external, so set them before you build.

**Which metabolites are the same compound?** Two units meet in the pool
only if their metabolites resolve to the same name. By default the
compartment is stripped, which resolves ``glc_D_e`` and ``glc_D[e]`` to the
same ``glc_D``, so the two spellings share one pool row in either mode.

What it cannot do is reconcile identifiers from *different namespaces*.
``glc__D_e`` and ``glc_D[e]`` are the same molecule under two schemes, and
no amount of suffix-stripping will tell you so. Left alone that failure is
silent and total: every unit trades through a pool of its own and nothing
cross-feeds, with no error anywhere. So the builder checks, and says so:

.. code-block:: text

   WARNING pool 'u' is shared by 2 units but not one metabolite in it is
   reached by more than one of them, so nothing can cross-feed. The usual
   cause is identifiers from different namespaces; pass metabolite_key to
   say how they correspond. Examples: {'A': ['glc__D_e'], 'B': ['glc_D[e]']}

``metabolite_key`` is the way out. Given a metabolite, return the name of
the thing it *is*:

.. code-block:: python

   import re

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
   :caption: ``A`` can reach ``b`` and export it but cannot finish; ``B``
             can finish but cannot start. Dashed reactions are the ones
             each draft is missing.
   :align: center

   digraph chain {
     rankdir=LR;
     bgcolor="transparent";
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#555555", arrowsize=0.7];

     node [shape=ellipse, style=filled, fillcolor="#e8f0fa",
           color="#3c6997", fontcolor="#1b3a57"];
     glc; a; b; c;

     node [shape=box, width=0.3, height=0.22, style=filled,
           fillcolor="#f7f7f7", color="#999999", fontcolor="#333333"];
     UP; R1; R2; BIO;

     glc -> UP -> a -> R1 -> b -> R2 -> c -> BIO;
     b -> EXP [label=" to the pool", fontcolor="#c07d2a"];
     EXP [shape=box, width=0.3, height=0.22, style=filled,
          fillcolor="#fdf0e3", color="#c07d2a", fontcolor="#6b4411"];
   }

Build the community. Each unit is backed by a *database* — everything it is
allowed to borrow — and its own draft supplies the bounds:

.. code-block:: python

   from spectra_cobra import build_community_model

   community = build_community_model(
       [draft_a, draft_b],
       organisms=["A", "B"],
       databases={"A": chain, "B": chain},
       pool_medium={"glc_e": (-10.0, 1000.0)},
   )

   community.database_reactions["A"]   # ('R2__A',)
   community.database_reactions["B"]   # ('UP__B', 'R1__B')

``pool_medium`` is what the environment offers, keyed by metabolite as your
own models spell it. Everything unnamed stays shut to uptake.

Now fill the gaps in both at once:

.. code-block:: python

   from spectra_cobra import gapfill_community

   result = gapfill_community(community, problem_type="minNetMILP")

   result.added        # {'A': ('R2',), 'B': ()}
   result.models()     # one untagged cobra model per unit

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
   * - already in the unit's draft
     - 0
     - it is not an addition
   * - only in the unit's database
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
exchange carries the flux of every unit drawing on it. Charging for it
therefore penalises a unit for having company, and the solver responds by
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
exchanges between them — the pool has one set of them however many units
there are, so the per-unit charge falls as the community grows. It is also
slower while it is wrong: 2,931 seconds against 854 at six organisms.

The mixed-integer formulation does not have the problem, because its
objective counts reactions rather than flux, so an exchange costs the same
whether one unit uses it or twenty. At two organisms it adds 8 against 10
filled separately, with the weight left at 1.

So the default follows the formulation: 1 for the MILP, 0 for the LP. Pass
``exchange_weight`` to override it either way. The published pipeline uses
1 throughout, and at the scale it reports on — a hundred and four
organisms, with the mixed-integer formulation — that is the right number.
It is small communities under the LP where it misleads.

Two more defaults worth knowing. Each unit's anchor reaction is **core**, so
every unit has to work — otherwise the cheapest community is one where half
the members are dead. And ``keep_draft=True`` keeps every reaction a unit
already had: the formulation returns the *smallest network* meeting the
requirement, which left to itself would throw away parts of the draft that
carry no flux under this particular medium. Gap-filling adds; it should not
quietly subtract.

Coupling, and what a unit is gated on
-------------------------------------

An absent unit must be absent in full. Not merely not growing — carrying no
flux at all. Otherwise a dead organism goes on running its metabolism and
feeding its neighbours for free.

That is what coupling does. For every reaction of a unit,

.. math::

   |v_i| \le c \cdot v_{\text{anchor}} + u

so a reaction can carry flux only in proportion to its unit's **anchor**.
The defaults, :math:`c = 1000` and :math:`u = 0.01`, are the published
values for a microbial community whose anchor is biomass.

Which reaction is the anchor is decided in this order, and never guessed
past it:

1. what you passed in ``biomass_reactions``,
2. the model's objective, if exactly one reaction is in it,
3. a single reaction named ``biomass*``,
4. otherwise it refuses, and tells you what the candidates looked like.

Refusing matters. A tissue is rarely growing, so its anchor is usually ATP
maintenance or a demand reaction with no telling name, and coupling ten
thousand reactions to the wrong one would throttle or free the whole tissue
without raising anything.

Whatever is chosen is reported:

.. code-block:: python

   print(community.coupling_summary())

.. code-block:: text

   unit    anchor                  how          capacity        cap
   tis1    DM_atp_c___tis1         given             10.0      1e+04
   tis2    DM_atp_c___tis2         given           0.0001       0.11

``capacity`` is the most that anchor can carry; ``cap`` is what the coupling
then permits every reaction in the unit. The second row is the failure to
watch for. An anchor running at 1e-4 with :math:`c = 1000` caps the whole
tissue at 0.11, which strangles it silently. The builder warns when it sees
this.

.. note::

   :math:`c` is not a number to be derived; it is a modelling choice about
   how much flux a unit may carry per unit of anchor. If your anchor runs at
   a different scale from a bacterial growth rate, raise :math:`c` to match
   it, or pass ``couple=False`` and do without.

A multi-tissue model
--------------------

Tissues need the ``"pooled"`` mode, for two reasons: they do not all meet in
the same place, and measured uptake rates attach to the transports.

Suppose three tissues. ``tis1`` and ``tis2`` are perfused by blood; ``tis3``
only ever meets ``tis2``:

.. graphviz::
   :caption: Two pools. ``Bl`` reaches the environment; the ``tis2_tis3``
             interface does not, so ``tis3`` lives entirely on what
             ``tis2`` passes it.
   :align: center

   digraph tissues {
     rankdir=LR;
     bgcolor="transparent";
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#555555", arrowsize=0.7,
           dir=both];

     node [shape=box, style="filled,rounded", fillcolor="#e8f0fa",
           color="#3c6997", fontcolor="#1b3a57"];
     tis1; tis2; tis3;

     node [shape=ellipse, style=filled, fillcolor="#fdf0e3",
           color="#c07d2a", fontcolor="#6b4411"];
     Bl [label="Bl"]; iface [label="tis2_tis3"];

     node [shape=point, width=0.05, color="#bbbbbb"];
     env;

     tis1 -> Bl; tis2 -> Bl;
     tis2 -> iface; tis3 -> iface;
     Bl -> env [style=dashed, color="#bbbbbb", label=" exchanges"];
   }

.. code-block:: python

   from spectra_cobra import build_multi_tissue_model

   body = build_multi_tissue_model(
       gem,
       tissues=["tis1", "tis2", "tis3"],
       pools={"Bl": ["tis1", "tis2"], "tis2_tis3": ["tis2", "tis3"]},
       environment=["Bl"],
       anchor_reactions={t: "DM_atp_c_" for t in ("tis1", "tis2", "tis3")},
       pool_medium=blood_medium,
       link_bounds=measured_rates,
       couple=False,
   )

Four things to notice.

``environment=["Bl"]`` is required, not optional. With more than one pool
the builder will not guess, because guessing "all of them" would give the
``tis2_tis3`` interface its own exchanges and let those two tissues bypass
the blood entirely.

``anchor_reactions`` names ATP maintenance. A tissue is not growing, so
there is no biomass reaction to find.

``link_bounds`` is where blood metabolomics goes. It bounds the transport
between a tissue and a pool — *this tissue may take up at most this much* —
which is a different statement from ``pool_medium``, which says what the
blood contains at all.

``couple=False`` here. With an anchor running at maintenance rates the
default :math:`c` would throttle every tissue; either raise :math:`c` to
suit the anchor or leave the coupling off. Build it, read
``coupling_summary()``, and decide with the numbers in front of you.

The reference network is read, never modified, and never copied: the same
model is handed over once per tissue, so a five-tissue human model costs no
more memory than a one-tissue one.

What that comes to, on a real reconstruction
--------------------------------------------

Four tissues on a consistent Recon3D, three on the blood and a fourth
reachable only through the third:

.. code-block:: text

   body: 36180 reactions, 20552 metabolites        (28s)
     pool Bl:        ['tis1', 'tis2', 'tis3'], 880 metabolites
     pool tis3_tis4: ['tis3', 'tis4'],         880 metabolites
     community exchanges: 880
     transports:         4400

4,400 is 880 × 3 for the blood plus 880 × 2 for the interface. The 880
exchanges are all on ``Bl``, because that is the only pool nominated.

The test that the topology is really doing something:

.. code-block:: text

   whole-body objective:       21.3061
   tis4 alone:                 21.3061
   tis4 with tis3 shut down:   0

``tis4`` touches no open pool, so everything it eats comes through the
interface from ``tis3``. Shut ``tis3`` down and ``tis4`` goes to exactly
zero. Had the interface been given exchanges of its own — which is what
defaulting ``environment`` to every pool would have done — ``tis4`` would
have carried on regardless and nothing would have looked wrong.

Extracting context-specific tissues
-----------------------------------

A multi-tissue model is a universal model like any other, so
:func:`~spectra_cobra.spectra_me` extracts from it directly. The core set is
the union of each tissue's core reactions, tagged:

.. code-block:: python

   from spectra_cobra import spectra_me

   core = [f"{rxn}__{tissue}" for tissue, rxns in cores.items() for rxn in rxns]
   extracted = spectra_me(body.model, core, tol=1e-4, problem_type="minNetLP")

   tissues = body.with_model(extracted).decompose()

On the four-tissue body above, with ATP maintenance and a maintenance
biomass reaction core in every tissue, ``minNetLP`` takes 36,180 reactions
down to 2,584 in about fourteen minutes, keeping all eight core reactions:

.. code-block:: text

   tis1 690   tis2 716   tis3 741   tis4 430     core kept: 8/8

``tis3`` comes back with the most boundary reactions of the four, which is
what you would hope: it is the only tissue sitting on two pools.

``decompose`` hands each tissue back as a model in its own right, untagged,
with an exchange reaction wherever it met a pool — a transport ``x[e] → x_u``
with the pool side dropped *is* an exchange, and it keeps the transport's
bounds, so the measured rates survive the round trip.

Community gap-filling at genome scale
-------------------------------------

The published case study gap-fills the hCom synthetic gut community against
CarveMe's universal reconstructions. Each organism's compartment holds an
entire universal — the Gram-positive, Gram-negative or common one, by Gram
stain — bounded by its own draft wherever the two agree. The organisms share
the extracellular compartment, constrained by the standard amino acid
complete (SAAC) medium.

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

Reference
---------

* :doc:`../functions/community` — the builder and the minimal microbiome
* :doc:`../functions/gapfilling` — ``gapfill_community`` and its siblings
* :doc:`minimal_microbiome` — reducing a community instead of filling it
