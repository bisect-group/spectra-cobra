Multi-tissue models
===================

A multi-tissue model is several copies of a metabolic network, one per
tissue, joined by the places they exchange metabolites: blood, an interstitial
space, the gut lumen. Each tissue keeps its own copy of every reaction and
every metabolite, so the same reaction running in liver and in muscle is two
separate reactions carrying separate flux. What connects them is a **pool** —
a compartment several tissues can reach — and what distinguishes them is
their bounds, which pools they reach, and which reactions they are required
to keep.

This page builds one from Recon3D and extracts context-specific tissues from
it. Nothing here assumes you have read anything else.

How the pieces are named
------------------------

Everything a tissue owns is tagged with that tissue's name. A reaction
``PGI`` becomes ``PGI__tis1`` in the first tissue and ``PGI__tis2`` in the
second; the external metabolite ``glc_D[e]`` becomes ``glc_D[e]__tis1`` and
``glc_D[e]__tis2``. The tag is a ``__`` suffix, which is why tissue names may
not contain ``__`` themselves.

A pool metabolite is **not** tagged, because it belongs to no tissue: blood
glucose is ``glc_D_Bl``, one metabolite that every tissue on the blood can
reach. Between a tissue's own external metabolite and the pool's sits a
**transport** reaction, one per tissue per metabolite per pool, named
``TR_<metabolite>_<pool>__<tissue>``. The pool reaches the outside world
through an **exchange**, ``EX_<metabolite>_<pool>``.

That arrangement — a private external compartment per tissue, and a transport
into each pool it reaches — is what lets one tissue belong to two pools at
once, and it is what gives you somewhere to attach a measured uptake rate.

.. graphviz::
   :caption: Three tissues and two pools. ``tis2`` is in both, and can be
             because it keeps its own copy of glucose and runs a separate
             transport into each pool. ``Bl`` reaches the environment;
             the ``tis2``/``tis3`` interface does not. Every ``TR_``
             belongs to exactly one tissue, which is what lets a measured
             rate attach to it.
   :align: center

   digraph pooled {
     rankdir=LR; bgcolor="transparent"; compound=true; nodesep=0.30;
     ranksep=0.55;
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#6b7280", arrowsize=0.7];

     subgraph cluster_a {
       label="tis1"; labeljust=l; fontname="Helvetica"; fontsize=9;
       fontcolor="#3c6997"; color="#c7d6e8"; style="rounded,filled";
       fillcolor="#f6f9fd"; margin=10;
       glcA [label="glc_t1", shape=ellipse, style=filled, fillcolor="#e8f0fa",
             color="#3c6997", fontcolor="#1b3a57"];
       trA1 [label="TR_t1_Bl", shape=box, style="filled,rounded", width=1.0,
             height=0.26, fillcolor="#ffffff", color="#3c6997",
             fontcolor="#1b3a57"];
       glcA -> trA1 [dir=both];
     }

     subgraph cluster_b {
       label="tis2  —  in both pools"; labeljust=l; fontname="Helvetica";
       fontsize=9; fontcolor="#b3411f"; color="#e2bfaf";
       style="rounded,filled"; fillcolor="#fdf6f3"; margin=10;
       glcB [label="glc_t2", shape=ellipse, style=filled, fillcolor="#e8f0fa",
             color="#3c6997", fontcolor="#1b3a57"];
       trB1 [label="TR_t2_Bl", shape=box, style="filled,rounded", width=1.0,
             height=0.26, fillcolor="#ffffff", color="#3c6997",
             fontcolor="#1b3a57"];
       trB2 [label="TR_t2_if", shape=box, style="filled,rounded", width=1.0,
             height=0.26, fillcolor="#ffffff", color="#3c6997",
             fontcolor="#1b3a57"];
       glcB -> trB1 [dir=both];
       glcB -> trB2 [dir=both];
     }

     subgraph cluster_c {
       label="tis3"; labeljust=l; fontname="Helvetica"; fontsize=9;
       fontcolor="#3c6997"; color="#c7d6e8"; style="rounded,filled";
       fillcolor="#f6f9fd"; margin=10;
       glcC [label="glc_t3", shape=ellipse, style=filled, fillcolor="#e8f0fa",
             color="#3c6997", fontcolor="#1b3a57"];
       trC2 [label="TR_t3_if", shape=box, style="filled,rounded", width=1.0,
             height=0.26, fillcolor="#ffffff", color="#3c6997",
             fontcolor="#1b3a57"];
       glcC -> trC2 [dir=both];
     }

     u1 [label="glc_Bl", shape=ellipse, style=filled, fillcolor="#fdf0e3",
         color="#c07d2a", fontcolor="#6b4411", penwidth=1.6];
     u2 [label="glc_if", shape=ellipse, style=filled, fillcolor="#fdf0e3",
         color="#c07d2a", fontcolor="#6b4411", penwidth=1.6];
     ex [label="EX_glc_Bl", shape=box, style="filled,rounded", width=1.0,
         height=0.26, fillcolor="#ffffff", color="#c07d2a",
         fontcolor="#6b4411"];
     env [label="environment", shape=plaintext, fontcolor="#9ca3af",
          fontsize=9];
     shut [label="no exchange —\nclosed to the outside", shape=plaintext,
           fontcolor="#b3411f", fontsize=8.5];

     trA1 -> u1 [dir=both, color="#c07d2a"];
     trB1 -> u1 [dir=both, color="#c07d2a"];
     trB2 -> u2 [dir=both, color="#c07d2a"];
     trC2 -> u2 [dir=both, color="#c07d2a"];
     u1 -> ex [dir=both, color="#c07d2a"];
     ex -> env [dir=both, style=dashed, color="#c9ced6"];
     u2 -> shut [dir=none, style=dotted, color="#d8b2a6"];
   }

The figure uses short labels to stay readable. In the model you build, those
nodes are ``glc_D[e]__tis1``, ``glc_D[e]__tis2``, ``glc_D[e]__tis3`` for the
three tissue copies, ``glc_D_Bl`` and ``glc_D_tis2_tis3`` for the two pools,
``TR_glc_D_Bl__tis1`` for a transport, and ``EX_glc_D_Bl`` for the exchange.

Building one
------------

Three tissues. ``tis1`` and ``tis2`` are perfused by blood; ``tis3`` only ever
meets ``tis2``.

.. graphviz::
   :caption: Two pools. ``Bl`` reaches the environment; the ``tis2_tis3``
             interface does not, so ``tis3`` lives on what ``tis2``
             passes it.
   :align: center

   digraph tissues {
     rankdir=LR; bgcolor="transparent"; nodesep=0.35; ranksep=0.55;
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#c07d2a", arrowsize=0.7,
           dir=both];

     node [shape=box, style="filled,rounded", fillcolor="#f6f9fd",
           color="#3c6997", fontcolor="#1b3a57", width=0.8, height=0.34];
     tis1; tis2; tis3;

     node [shape=ellipse, style=filled, fillcolor="#fdf0e3", color="#c07d2a",
           fontcolor="#6b4411", penwidth=1.6];
     Bl    [label="Bl"];
     iface [label="tis2_tis3"];

     ex  [label="EX_*_Bl", shape=box, style="filled,rounded", width=0.8,
          height=0.26, fillcolor="#ffffff", color="#c07d2a",
          fontcolor="#6b4411"];
     env [label="environment", shape=plaintext, style="", fillcolor="transparent",
          fontcolor="#9ca3af", fontsize=9];

     tis1 -> Bl;
     tis2 -> Bl;
     tis2 -> iface;
     tis3 -> iface;
     Bl -> ex;
     ex -> env [style=dashed, color="#c9ced6"];

     shut [label="no exchanges —\nclosed to the outside", shape=plaintext,
           style="", fillcolor="transparent", fontcolor="#b3411f", fontsize=8.5];
     iface -> shut [dir=none, style=dotted, color="#d8b2a6"];
     { rank=same; ex; shut; }
   }

The reference network is Recon3D, 10,600 reactions over 5,835 metabolites.
``Recon3DModel_301.mat`` is in the COBRA Toolbox's open model collection, and
every number on this page was measured against that exact file:

   https://github.com/opencobra/COBRA.models/blob/master/mat/Recon3DModel_301.mat

:doc:`recon3d` covers working with it at genome scale.

.. code-block:: python

   from cobra.io import load_matlab_model
   from spectra_cobra import build_multi_tissue_model

   gem = load_matlab_model("Recon3DModel_301.mat")

   # What the blood carries. Recon3D ships with all of its exchanges open,
   # so this says the blood offers everything the reconstruction knows how
   # to take up. Narrow it to the metabolites your blood actually carries
   # when you have that measurement -- but narrow it deliberately, because
   # a human tissue starves on a medium that an E. coli model would be
   # perfectly happy with.
   blood = {next(iter(r.metabolites)).id: (-1000.0, 1000.0)
            for r in gem.reactions
            if r.id.startswith("EX_") and len(r.metabolites) == 1}

   body = build_multi_tissue_model(
       gem,
       tissues=["tis1", "tis2", "tis3"],
       pools={"Bl": ["tis1", "tis2"], "tis2_tis3": ["tis2", "tis3"]},
       environment=["Bl"],
       anchor_reactions={t: "biomass_maintenance"
                         for t in ("tis1", "tis2", "tis3")},
       pool_medium=blood,
       couple=False,
   )

   print(f"{len(body.model.reactions)} reactions, "
         f"{len(body.model.metabolites)} metabolites")
   for name, units in body.pools.items():
       print(f"  pool {name}: {list(units)}, "
             f"{len(body.pool_metabolites[name])} metabolites")
   print(f"  community exchanges: {len(body.community_exchanges)}")
   transports = sum(r.id.startswith("TR_") for r in body.model.reactions)
   print(f"  transports: {transports}")

.. code-block:: text

   35015 reactions, 20663 metabolites
     pool Bl: ['tis1', 'tis2'], 1579 metabolites
     pool tis2_tis3: ['tis2', 'tis3'], 1579 metabolites
     community exchanges: 1579
     transports: 6316

6,316 is 1,579 × 2 for the blood plus 1,579 × 2 for the interface. The 1,579
exchanges are all on ``Bl``, because that is the only pool nominated.

Four things to notice.

``environment=["Bl"]`` is required, not optional. With more than one pool the
builder will not guess, because guessing "all of them" would give the
``tis2_tis3`` interface its own exchanges and let those two tissues feed
themselves from outside, bypassing the blood entirely.

``anchor_reactions`` names ``biomass_maintenance``. A tissue is not growing,
so there is no biomass reaction to find and the anchor has to be said rather
than guessed. Recon3D's maintenance biomass is the natural choice: it is what
the tissue must keep doing to stay alive.

``couple=False`` here. Coupling ties a tissue's reactions to its anchor, so
that nothing runs in a tissue that is not maintaining itself. With an anchor
running at maintenance rates the default coupling factor would throttle every
tissue, so either raise it to suit the anchor or leave the coupling off.
Build it, read ``coupling_summary()``, and decide with the numbers in front
of you.

``pool_medium`` is what the pool contains at all. It is a different statement
from ``link_bounds``, which bounds a single tissue's transport into a pool —
*this* tissue may take up at most this much. Blood metabolomics goes in the
first; measured per-tissue uptake rates go in the second.

Tissues that are not copies of each other
-----------------------------------------

Above, every tissue starts as a copy of one reconstruction, which is the
usual case: what makes liver differ from muscle is its bounds and its core
reactions, not a different reaction database. When the tissues really are
different networks — already context-specific models, say, or models from
different sources — hand over one per tissue instead, in the same order:

.. code-block:: python

   body = build_multi_tissue_model(
       [liver_model, muscle_model, adipose_model],
       tissues=["liver", "muscle", "adipose"],
       pools={"Bl": ["liver", "muscle", "adipose"]},
       environment=["Bl"],
       anchor_reactions={"liver": "biomass_maintenance",
                         "muscle": "biomass_maintenance",
                         "adipose": "biomass_maintenance"},
   )

Everything else is unchanged. Giving one model costs no more memory than a
one-tissue model would: the builder reads its inputs and never modifies them,
so the same object is handed over once per tissue rather than copied.

Does the topology actually bite?
--------------------------------

``tis3`` touches no open pool. Everything it gets should arrive through the
interface from ``tis2``, so shutting ``tis2`` down should cost ``tis3``
everything. Testing that is worth the two lines it takes, because a topology
that is not doing what you think looks exactly like one that is:

.. code-block:: python

   body.model.objective = body.model.reactions.get_by_id(
       "biomass_maintenance__tis3")
   print(f"tis3 alone:               {body.model.slim_optimize():.4f}")

   with body.model:
       for rxn_id in body.reactions_of["tis2"]:
           body.model.reactions.get_by_id(rxn_id).bounds = (0.0, 0.0)
       print(f"tis3 with tis2 shut down: "
             f"{body.model.slim_optimize(error_value=0.0):.4f}")

.. code-block:: text

   tis3 alone:               755.0032
   tis3 with tis2 shut down: 150.8679

Shutting ``tis2`` down costs ``tis3`` four fifths of its maintenance capacity,
which is the interface doing its job. But it does not cost everything, and the
reason is worth knowing.

.. warning::

   **A human reconstruction ships with taps no pool controls.** Recon3D has
   247 boundary reactions that are not named ``EX_``: 151 demands (``DM_``)
   and 96 sinks (``sink_``). Every one of the 96 sinks has a negative lower
   bound, meaning it supplies a metabolite from nothing; none of the demands
   does. They sit *inside* the cell, so tagging gives every tissue its own —
   ``sink_crvnc[c]__tis3`` is an open tap in ``tis3`` that neither blood nor
   the interface governs.

   That floor of 150.87 is what those sinks provide. Close them and the same
   body behaves as the figure claims:

   .. list-table::
      :header-rows: 1
      :widths: 34 16 16 16 18

      * - sinks
        - ``tis1``
        - ``tis2``
        - ``tis3``
        - ``tis3``, ``tis2`` shut
      * - as Recon3D ships them
        - 755.0032
        - 755.0032
        - 755.0032
        - 150.8679
      * - 96 supplying sinks closed
        - 684.3324
        - **684.5431**
        - 684.3324
        - **0.0000**

   The second row is the more interesting one. With the sinks open all three
   tissues read an identical 755.0032: whether a tissue sits on one pool or
   two makes no difference, because the sinks dominate. Closed, ``tis2``
   separates out — 684.5431 against 684.3324 — and that small margin is
   exactly the advantage of being the only tissue on both pools. The sinks do
   not merely soften the result; they erase the structural signal the model
   was built to express, while leaving a number that looks perfectly
   reasonable.

   Closing them is one line, and whether to is a modelling decision rather
   than a default:

   .. code-block:: python

      for rxn in gem.reactions:
          if rxn.id.startswith(("DM_", "SK_", "sink_")) and rxn.lower_bound < 0:
              rxn.lower_bound = 0.0

Where each tissue's inputs come from
-------------------------------------

A multi-tissue model has four inputs, and in real work each comes from a
different measurement.

**The pool medium** is what the blood carries, from blood metabolomics. It
sets what is available at all. Above it is left wide open, which keeps the
page about topology — but a wide-open blood is a strong assumption, and
narrowing it is usually the first thing to do with real data.

**The link bounds** are what each tissue may draw from a pool, from measured
per-tissue uptake and secretion rates. ``link_bounds`` attaches them to the
transports, which is why the builder makes one transport per tissue rather
than letting tissues share.

**The anchor** is what each tissue must keep doing. Maintenance biomass, or
ATP maintenance, or a demand reaction for a product the tissue is known to
export.

**The core reactions** are what each tissue is known to be running, and this
is where transcriptomics goes. The usual route is to call a gene present or
absent from expression, map genes to reactions through the model's GPR rules,
and take the reactions whose rules are satisfied as that tissue's core. Those
reaction identifiers are then tagged with the tissue's name — ``PGI__tis1`` —
and the tagged sets from every tissue are concatenated into one core list for
the whole body.

.. note::

   The example below uses each tissue's maintenance reaction as its core,
   chosen by hand. That keeps the page about the mechanics. In real work the
   core comes from expression data as described above, and it is a great deal
   larger: a few hundred to a few thousand reactions per tissue rather than
   one. The extraction call does not change — only the contents of ``core``.

Extracting context-specific tissues
-----------------------------------

A multi-tissue model is a universal model like any other, so
:func:`~spectra_cobra.spectra_me` extracts from it directly. The core set is
the union of each tissue's core reactions, tagged:

.. code-block:: python

   from spectra_cobra import spectra_me, check_extraction

   # One reaction per tissue here; from transcriptomics it would be many.
   core = [f"biomass_maintenance__{t}" for t in ("tis1", "tis2", "tis3")]
   extracted = spectra_me(body.model, core, tol=1e-4,
                          problem_type="minNetLP", seed=0)

   tissues = body.with_model(extracted).decompose()

   print(f"{len(extracted.reactions)} reactions")
   print("  " + "   ".join(f"{name} {len(model.reactions)}"
                           for name, model in tissues.items()))
   print(f"  {check_extraction(extracted, core, tol=1e-4)}")

.. code-block:: text

   357 reactions
     tis1 97   tis2 115   tis3 125
     <ExtractionReport 357 reactions, 0 core missing, 0 blocked (0 of them core), valid=True>

``decompose`` splits the extracted body back into one model per tissue, each
untagged and standing on its own, with an exchange reaction wherever it met a
pool — a transport with the pool side dropped *is* an exchange, and it keeps
the transport's bounds, so measured rates survive the round trip.

Always run :func:`~spectra_cobra.check_extraction`. A core reaction being
*present* in the result is not the same as its being able to carry flux, and
only the check distinguishes them.

.. warning::

   ``minNetLP`` infers which reactions it kept from which ones carry flux in
   its solution. That inference is sound only while the fluxes mean
   something. The objective pushes every reaction as low as the constraints
   allow, and on a large body a share of them can end up carrying less than
   the solver's own feasibility tolerance — numbers it cannot distinguish
   from zero. Reactions kept on that evidence may leave the extracted model
   unable to carry flux at all, core reactions included.

   SPECTRA warns when it sees the condition:

   .. code-block:: text

      WARNING N of M reactions in this solution carry less flux than the
      solver can resolve (model.tolerance=1e-09), so whether they belong in
      the model is not something the solution settles. The extracted model
      may be unable to carry flux at all, core reactions included.

   The remedy is a larger ``tol``, which lifts the whole solution clear of
   the solver's resolution. No inclusion cutoff helps, because the
   information was never in the solution. How much larger is not something
   that can be said in general — raise it until the warning stops, and
   confirm with :func:`~spectra_cobra.check_extraction`.

   The three-tissue body above does not trigger it at ``tol=1e-4``: all 357
   reactions carry flux the solver can resolve, and the result verifies.
   Larger bodies, more tissues and larger cores all make it more likely.

Other formulations, and what comes after
----------------------------------------

Nothing in that call is specific to a multi-tissue body. ``spectra_me``
takes a universal model and a core set, and a tagged body is simply a larger
universal model: run the same call on a single untagged GEM with one
tissue's core reactions and you have the ordinary single-model
context-specific extraction that :doc:`extraction` describes. The tags are
only what lets one solve cover three tissues while keeping each tissue's
flux its own.

``problem_type`` is equally free. ``minNetLP`` is used above because it is
cheap — one LP, 29 s on this 35,015-reaction body — but it is one of four,
and they ask genuinely different questions rather than rounding one answer
four ways:

.. list-table::
   :header-rows: 1
   :widths: 24 46 15 15

   * - ``problem_type``
     - Objective
     - Weights
     - Solves
   * - ``"minNetLP"``
     - minimise weighted total absolute flux
     - non-negative
     - 1 LP
   * - ``"minNetMILP"``
     - minimise weighted reaction count
     - non-negative
     - 1 MILP
   * - ``"growthOptim"``
     - maximise objective flux minus weighted flux
     - non-negative
     - 1 LP
   * - ``"tradeOff"``
     - maximise net weight: positive in, negative out
     - any real
     - 1 MILP

:doc:`formulations` sets out when each is the right question. The two
mixed-integer ones answer a sharper version of it — *fewest reactions*
rather than *least flux* — and on a body this size they are worth giving an
explicit ``time_limit`` and reading as "best found so far" rather than as a
proven optimum.

Extraction also need not be the last step. A model built from expression
evidence can still be unable to do things the tissue demonstrably does,
because the evidence never covered them. Metabolic tasks state those
capabilities directly — regenerate ATP from glucose aerobically, synthesise
this amino acid from that precursor — and
:func:`~spectra_cobra.check_tasks` reports which ones a tissue fails.
:func:`~spectra_cobra.gapfill_for_tasks` then adds back the fewest reactions
from a universal model that make the failing ones pass, with weights
carrying the same evidence that drove the extraction. It applies on either
side of ``decompose``: per tissue, with the untagged GEM as the donor, or to
the whole body at once, where the candidates include the transports — so a
capability one tissue lacks can be restored by letting it draw the
metabolite from a pool rather than by giving it the pathway.
:doc:`gapfilling_tasks` works through both the task format and the solve.

References
----------

* Brunk, E., Sahoo, S., Zielinski, D. C., et al. (2018). Recon3D enables a
  three-dimensional view of gene variation in human metabolism. *Nature
  Biotechnology*, 36(3), 272–281. https://doi.org/10.1038/nbt.4072

See also
--------

* :doc:`../functions/community` — the builder's full signature
* :doc:`recon3d` — working at genome scale, and where to get the model
* :doc:`extraction` — what :func:`~spectra_cobra.spectra_me` does
* :doc:`formulations` — choosing a ``problem_type``
* :doc:`gapfilling_tasks` — holding a model to what it must be able to do
