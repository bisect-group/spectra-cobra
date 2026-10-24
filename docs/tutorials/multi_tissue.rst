Multi-tissue models
===================

Tissues of a body are built exactly as the organisms of a community are: one
tagged copy of the reference network per tissue, meeting in pools. Everything
in :doc:`community` about pools, tagging, exchangeable metabolites and
coupling applies here unchanged, and this page assumes it.

Three things differ in practice, and all three come from the fact that a
tissue is not an organism:

* **Every tissue is a copy of the same network.** What distinguishes them is
  their bounds, their core reactions and which pool they can reach, so
  :func:`~spectra_cobra.build_multi_tissue_model` takes one model and a list
  of names rather than a list of models.
* **A tissue is rarely growing.** There is no biomass reaction to find, so
  the anchor that coupling gates on has to be named — usually maintenance or
  a demand reaction.
* **The topology is rarely one pool.** Blood reaches some tissues; others
  meet only each other. That is what ``pools`` and ``environment`` are for.

A multi-tissue model
--------------------

Tissues need the ``"pooled"`` mode, for two reasons: they do not all meet in
the same place, and measured uptake rates attach to the transports.

What ``"pooled"`` builds is a private external compartment per tissue and a
transport reaction from it into each pool that tissue reaches. That is what
lets one tissue belong to more than one pool at once:

.. graphviz::
   :caption: ``mode="pooled"``, with two pools. ``B`` is in both, and can be
             because it keeps its own ``glc[e]`` and runs a separate
             transport into each. ``u1`` reaches the environment; ``u2`` is an
             interface between ``B`` and ``C`` and stays closed. Each ``TR_``
             belongs to one tissue, which is what ``link_bounds`` attaches a
             measured rate to.
   :align: center

   digraph pooled {
     rankdir=LR; bgcolor="transparent"; compound=true; nodesep=0.30;
     ranksep=0.55;
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#6b7280", arrowsize=0.7];

     subgraph cluster_a {
       label="tissue A"; labeljust=l; fontname="Helvetica"; fontsize=9;
       fontcolor="#3c6997"; color="#c7d6e8"; style="rounded,filled";
       fillcolor="#f6f9fd"; margin=10;
       glcA [label="glc[e]", shape=ellipse, style=filled, fillcolor="#e8f0fa",
             color="#3c6997", fontcolor="#1b3a57"];
       trA1 [label="TR_glc_u1", shape=box, style="filled,rounded", width=1.0,
             height=0.26, fillcolor="#ffffff", color="#3c6997",
             fontcolor="#1b3a57"];
       glcA -> trA1 [dir=both];
     }

     subgraph cluster_b {
       label="tissue B  —  in both pools"; labeljust=l; fontname="Helvetica";
       fontsize=9; fontcolor="#b3411f"; color="#e2bfaf";
       style="rounded,filled"; fillcolor="#fdf6f3"; margin=10;
       glcB [label="glc[e]", shape=ellipse, style=filled, fillcolor="#e8f0fa",
             color="#3c6997", fontcolor="#1b3a57"];
       trB1 [label="TR_glc_u1", shape=box, style="filled,rounded", width=1.0,
             height=0.26, fillcolor="#ffffff", color="#3c6997",
             fontcolor="#1b3a57"];
       trB2 [label="TR_glc_u2", shape=box, style="filled,rounded", width=1.0,
             height=0.26, fillcolor="#ffffff", color="#3c6997",
             fontcolor="#1b3a57"];
       glcB -> trB1 [dir=both];
       glcB -> trB2 [dir=both];
     }

     subgraph cluster_c {
       label="tissue C"; labeljust=l; fontname="Helvetica"; fontsize=9;
       fontcolor="#3c6997"; color="#c7d6e8"; style="rounded,filled";
       fillcolor="#f6f9fd"; margin=10;
       glcC [label="glc[e]", shape=ellipse, style=filled, fillcolor="#e8f0fa",
             color="#3c6997", fontcolor="#1b3a57"];
       trC2 [label="TR_glc_u2", shape=box, style="filled,rounded", width=1.0,
             height=0.26, fillcolor="#ffffff", color="#3c6997",
             fontcolor="#1b3a57"];
       glcC -> trC2 [dir=both];
     }

     u1 [label="glc_u1", shape=ellipse, style=filled, fillcolor="#fdf0e3",
         color="#c07d2a", fontcolor="#6b4411", penwidth=1.6];
     u2 [label="glc_u2", shape=ellipse, style=filled, fillcolor="#fdf0e3",
         color="#c07d2a", fontcolor="#6b4411", penwidth=1.6];
     ex [label="EX_glc_u1", shape=box, style="filled,rounded", width=1.0,
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

A community is normally built the other way, with ``mode="shared"``: one
external compartment for everyone, no transports, and therefore exactly one
pool. :doc:`community` has that figure. The trade is size — one transport
per tissue per metabolite is affordable across five tissues and is not
across a hundred organisms.

Suppose three tissues. ``tis1`` and ``tis2`` are perfused by blood; ``tis3``
only ever meets ``tis2``:

.. graphviz::
   :caption: Two pools. ``Bl`` reaches the environment; the ``tis2_tis3``
             interface does not, so ``tis3`` lives entirely on what
             ``tis2`` passes it.
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

.. code-block:: python

   from cobra.io import load_model
   from spectra_cobra import build_multi_tissue_model

   # The reference network, and the two sets of bounds the body is built with.
   # Substitute your own: ``pool_medium`` is what the blood offers, and
   # ``link_bounds`` what may cross between a tissue and a pool.
   gem = load_model("iJO1366")   # or Recon3D, converted -- see :doc:`recon3d`
   blood_medium = {"glc__D_e": (-1.0, 1000.0), "o2_e": (-1000.0, 1000.0)}
   measured_rates = {"glc__D_e": (-0.5, 1000.0)}

   body = build_multi_tissue_model(
       gem,
       tissues=["tis1", "tis2", "tis3"],
       pools={"Bl": ["tis1", "tis2"], "tis2_tis3": ["tis2", "tis3"]},
       environment=["Bl"],
       anchor_reactions={t: "ATPM" for t in ("tis1", "tis2", "tis3")},
       pool_medium=blood_medium,
       link_bounds=measured_rates,
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

   8397 reactions, 6063 metabolites
     pool Bl: ['tis1', 'tis2'], 324 metabolites
     pool tis2_tis3: ['tis2', 'tis3'], 324 metabolites
     community exchanges: 324
     transports: 1296

*E. coli* standing in for a tissue is biologically odd, but the structure is
the point and it loads anywhere. 1296 is 324 x 2 for the blood plus 324 x 2
for the interface; the 324 exchanges are all on ``Bl``, because that is the
only pool nominated.

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

The test that the topology is really doing something, on the *E. coli* body
built above — where ``tis3`` is the tissue reachable only through ``tis2``:

.. code-block:: python

   body.model.objective = body.model.reactions.get_by_id("ATPM__tis3")
   print(f"tis3 alone:               {body.model.slim_optimize():.4f}")

   with body.model:
       for rxn_id in body.reactions_of["tis2"]:
           body.model.reactions.get_by_id(rxn_id).bounds = (0.0, 0.0)
       print(f"tis3 with tis2 shut down: "
             f"{body.model.slim_optimize(error_value=0.0):.4f}")

.. code-block:: text

   tis3 alone:               17.2000
   tis3 with tis2 shut down: 0.0000

``tis3`` touches no open pool, so everything it eats comes through the
interface from ``tis2``. Shut ``tis2`` down and ``tis3`` goes to exactly
zero. Had the interface been given exchanges of its own — which is what
defaulting ``environment`` to every pool would have done — ``tis3`` would
have carried on regardless and nothing would have looked wrong.

The Recon3D figures above show the same structure at genome scale, where
``tis4`` stands in the position ``tis3`` holds here.

Extracting context-specific tissues
-----------------------------------

A multi-tissue model is a universal model like any other, so
:func:`~spectra_cobra.spectra_me` extracts from it directly. The core set is
the union of each tissue's core reactions, tagged:

.. code-block:: python

   from spectra_cobra import spectra_me

   # ATP maintenance in every tissue. On a human reconstruction a maintenance
   # biomass reaction would sit alongside it in each tissue's core.
   core = [f"ATPM__{tissue}" for tissue in ("tis1", "tis2", "tis3")]
   extracted = spectra_me(body.model, core, tol=1e-4, problem_type="minNetLP",
                          seed=0)

   tissues = body.with_model(extracted).decompose()

   print(f"{len(extracted.reactions)} reactions")
   print("  " + "   ".join(f"{name} {len(model.reactions)}"
                           for name, model in tissues.items()))
   kept = {rxn.id for rxn in extracted.reactions}
   print(f"  core kept: {sum(c in kept for c in core)}/{len(core)}")

.. code-block:: text

   113 reactions
     tis1 35   tis2 39   tis3 35
     core kept: 3/3

``decompose`` splits the extracted body back into one model per tissue. The
tissues differ because they reach different pools: ``tis2`` is on both the
blood and the interface and carries the most.

On the four-tissue body above, with ATP maintenance and a maintenance
biomass reaction core in every tissue:

.. list-table::
   :header-rows: 1
   :widths: 22 10 16 14 38

   * - Formulation
     - ``tol``
     - Reactions
     - Time
     - :func:`~spectra_cobra.check_extraction`
   * - ``minNetLP``
     - 1e-4
     - 2,596
     - 881s
     - **fails**: 2,292 blocked, 5 of them core
   * - ``minNetLP``
     - 1e-3
     - 894
     - 948s
     - passes
   * - ``minNetMILP``
     - 1e-4
     - **420**
     - 4,140s
     - passes

.. warning::

   ``minNetLP`` needs a larger ``tol`` here than a single organism does,
   and the usual 1e-4 is not enough. At that value it returned a model
   holding all eight core reactions and unable to carry flux through
   seven of them -- and only :func:`~spectra_cobra.check_extraction` says
   so, because a core reaction being *present* is not the same as its
   working.

   The reason is that an LP infers which reactions it kept from which
   ones carry flux, and a forced core flux of 1e-4 sits only five orders
   of magnitude above a solver tolerance of 1e-9. On a body this size the
   objective drives half the solution into that gap: **645 of the 1,269
   reactions it kept carried less flux than the solver can resolve**, and
   1,267 of them turned out to be blocked. Raising ``tol`` to 1e-3 lifts
   the whole solution clear -- nothing below tolerance, and the result
   verifies.

   No inclusion cutoff fixes this, because the information was never in
   the solution; sweeping it from 1e-6 down to zero changes nothing.
   SPECTRA warns when it sees the condition:

   .. code-block:: text

      WARNING 645 of 1269 reactions in this solution carry less flux than
      the solver can resolve (model.tolerance=1e-09) ...

   ``minNetMILP`` is immune, because its membership comes from its
   binaries rather than from flux, and it verifies at ``tol=1e-4``.

So the choice at this scale is between a verified 894 reactions in about
sixteen minutes and a verified 420 in about seventy. Note also that the
mixed-integer run did not *finish*: it stopped at its one-hour limit with
the usual warning that the answer may not be minimal, and that unproven
incumbent is still less than half the size of the LP's.

.. code-block:: text

   tis1 244   tis2 246   tis3 270   tis4 127     core kept: 8/8

``tis3`` comes back with the most boundary reactions of the four, which is
what you would hope: it is the only tissue sitting on two pools.

``decompose`` hands each tissue back as a model in its own right, untagged,
with an exchange reaction wherever it met a pool — a transport ``x[e] → x_u``
with the pool side dropped *is* an exchange, and it keeps the transport's
bounds, so the measured rates survive the round trip.

References
----------

* Brunk, E., Sahoo, S., Zielinski, D. C., et al. (2018). Recon3D enables a
  three-dimensional view of gene variation in human metabolism. *Nature
  Biotechnology*, 36(3), 272–281. https://doi.org/10.1038/nbt.4072

Reference
---------

* :doc:`community` — the construction these tissues are built with
* :doc:`../functions/community` — the builder itself
* :doc:`recon3d` — working at genome scale
