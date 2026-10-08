Topology or stoichiometry
=========================

Every routine takes ``consistency_type``, which chooses what "can carry flux"
means:

``"stoichiometry"`` (the default)
   The usual steady state, :math:`S v = 0`. Everything produced is consumed.

``"topology"``
   An accumulation condition, :math:`S v \ge 0`. Metabolites may pile up but
   may not appear from nothing.

Topology is the weaker requirement, so it never finds *fewer* consistent
reactions. It is what topology-based gap filling uses, where the question is
whether a route to a metabolite exists at all rather than whether the whole
network balances.

A network the steady state rejects
----------------------------------

This is figure 2a of the Meneco article, and the toy model SPECTRA's own
experiment uses. Note the coefficient ``n``:

.. code-block:: python

   from cobra import Model, Reaction
   from spectra_cobra import consistent_reaction_ids


   def topology_toy_model(n):
       model = Model(f"topology n={n}")
       ids = [f"r{i}" for i in range(1, 7)]
       model.add_reactions([Reaction(i) for i in ids])
       formulas = ["--> S", f"S --> b + {n} a", "b + a --> d",
                   "a --> c", "d + c --> T1", "T1 -->"]
       for rxn_id, formula in zip(ids, formulas):
           rxn = model.reactions.get_by_id(rxn_id)
           rxn.reaction = formula
           rxn.bounds = (0.0, 1000.0)
       return model


   for n in (1, 2, 3):
       model = topology_toy_model(n)
       steady, _ = consistent_reaction_ids(model, tol=1e-4,
                                           consistency_type="stoichiometry")
       accumulating, _ = consistent_reaction_ids(model, tol=1e-4,
                                                 consistency_type="topology")
       print(f"n={n}: stoichiometry {len(steady)}/6, topology "
             f"{len(accumulating)}/6")

.. graphviz::
   :caption: ``topology_toy_model(n)``. Boxes are reactions, ovals are
             metabolites. Everything is irreversible, and only ``r2``
             carries a coefficient other than one.
   :align: center

   digraph topology_toy {
     rankdir=LR; bgcolor="transparent"; nodesep=0.24; ranksep=0.40;
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#6b7280", arrowsize=0.65];

     node [shape=ellipse, style=filled, fillcolor="#e8f0fa", color="#3c6997",
           fontcolor="#1b3a57", width=0.34, height=0.30];
     S; a; b; c; d; T1;

     node [shape=box, style="filled,rounded", width=0.44, height=0.26,
           fillcolor="#ffffff", color="#3c6997", fontcolor="#1b3a57"];
     r1; r2; r3; r4; r5; r6;

     node [shape=plaintext, style="", fillcolor="transparent",
           fontcolor="#9ca3af", fontsize=9];
     src [label="medium"]; sink [label="out"];

     src -> r1 [style=dashed, color="#c9ced6"];
     r1 -> S;
     S  -> r2;
     r2 -> b;
     r2 -> a [label=" n ", fontcolor="#b3411f", color="#b3411f", penwidth=1.5];
     b  -> r3;
     a  -> r3;
     r3 -> d;
     a  -> r4;
     r4 -> c;
     d  -> r5;
     c  -> r5;
     r5 -> T1;
     T1 -> r6;
     r6 -> sink [style=dashed, color="#c9ced6"];
   }

.. code-block:: text

   n=1: stoichiometry 0/6, topology 6/6
   n=2: stoichiometry 6/6, topology 6/6
   n=3: stoichiometry 0/6, topology 6/6

Not a rounding artefact. Chase the mass balances through and they force

.. math::

   v_2 \cdot (n - 2) = 0

so at steady state the entire network is dead unless ``n`` is exactly 2.
Relax the balances to inequalities and every ``n`` works, because the surplus
``a`` is allowed to accumulate instead of having to be consumed exactly.

What topology will and will not rescue
--------------------------------------

The weaker condition is not a free pass. ``S v >= 0`` lets a metabolite pile
up; it does not let one appear from nothing. These two networks differ only in
the direction of ``r3``:

.. code-block:: python

   # C is produced by r3, consumed by nothing.
   produced = ["--> A", "A --> B", "B --> C", "B --> D", "D -->",
               "A --> E", "E --> F", "F --> D"]

   # C is consumed by r3, produced by nothing.
   consumed = ["--> A", "A --> B", "C --> B", "B --> D", "D -->",
               "A --> E", "E --> F", "F --> D"]

Every reaction here has one substrate and one product, so the metabolites
alone are enough to draw it. Only the arrow on ``r3`` differs:

.. graphviz::
   :caption: The only difference is which way ``r3`` points. Produced and
             never consumed, ``C`` is a dead end that accumulation
             rescues; consumed and never produced, it is a gap that
             nothing can.
   :align: center

   digraph cc_toys {
     bgcolor="transparent"; rankdir=LR; nodesep=0.22; ranksep=0.38;
     node [fontname="Helvetica", fontsize=10];
     edge [fontname="Helvetica", fontsize=9, color="#6b7280", arrowsize=0.65];

     subgraph cluster_produced {
       label="C produced, never consumed  —  rescued by accumulation";
       labeljust=l; fontname="Helvetica"; fontsize=9; fontcolor="#3c6997";
       color="#c7d6e8"; style="rounded,filled"; fillcolor="#f6f9fd"; margin=10;
       node [shape=plaintext, style="", fillcolor="transparent",
             fontcolor="#9ca3af", fontsize=8.5];
       in1 [label="medium"]; out1 [label="out"];
       node [shape=ellipse, style=filled, fillcolor="#e8f0fa", color="#3c6997",
             fontcolor="#1b3a57", width=0.36, height=0.32];
       A1 [label="A"]; B1 [label="B"]; D1 [label="D"];
       E1 [label="E"]; F1 [label="F"];
       C1 [label="C", fillcolor="#fdf2ef", color="#b3411f",
           fontcolor="#b3411f", penwidth=1.6];
       in1 -> A1 [label=" r1", style=dashed, color="#c9ced6"];
       A1 -> B1 [label=" r2"];
       B1 -> C1 [label=" r3", color="#b3411f", fontcolor="#b3411f",
                 penwidth=1.6];
       B1 -> D1 [label=" r4"];
       D1 -> out1 [label=" r5", style=dashed, color="#c9ced6"];
       A1 -> E1 [label=" r6"];
       E1 -> F1 [label=" r7"];
       F1 -> D1 [label=" r8"];
     }

     subgraph cluster_consumed {
       label="C consumed, never produced  —  a gap nothing rescues";
       labeljust=l; fontname="Helvetica"; fontsize=9; fontcolor="#3c6997";
       color="#c7d6e8"; style="rounded,filled"; fillcolor="#f6f9fd"; margin=10;
       node [shape=plaintext, style="", fillcolor="transparent",
             fontcolor="#9ca3af", fontsize=8.5];
       in2 [label="medium"]; out2 [label="out"];
       node [shape=ellipse, style=filled, fillcolor="#e8f0fa", color="#3c6997",
             fontcolor="#1b3a57", width=0.36, height=0.32];
       A2 [label="A"]; B2 [label="B"]; D2 [label="D"];
       E2 [label="E"]; F2 [label="F"];
       C2 [label="C", fillcolor="#fdf2ef", color="#b3411f",
           fontcolor="#b3411f", penwidth=1.6];
       in2 -> A2 [label=" r1", style=dashed, color="#c9ced6"];
       A2 -> B2 [label=" r2"];
       C2 -> B2 [label=" r3", color="#b3411f", fontcolor="#b3411f",
                 penwidth=1.6];
       B2 -> D2 [label=" r4"];
       D2 -> out2 [label=" r5", style=dashed, color="#c9ced6"];
       A2 -> E2 [label=" r6"];
       E2 -> F2 [label=" r7"];
       F2 -> D2 [label=" r8"];
     }
   }

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - model
     - stoichiometry
     - topology
   * - ``C`` produced, never consumed
     - 7 of 8 (``r3`` blocked)
     - **8 of 8**
   * - ``C`` consumed, never produced
     - 7 of 8 (``r3`` blocked)
     - **7 of 8** (still blocked)

The first is rescued: letting ``C`` accumulate is exactly what it needed. The
second is not, and cannot be. For that metabolite the balance reads
:math:`-v_3 \ge 0`, which with :math:`v_3 \ge 0` still forces
:math:`v_3 = 0`.

So if topology mode does not unblock a reaction you expected it to, check
whether the gap is a missing *consumer* (which it fixes) or a missing
*producer* (which it does not — you need to add the producing reaction).

Extraction in topology mode
---------------------------

``consistency_type`` works the same way on the extraction routines:

.. code-block:: python

   from spectra_cobra import spectra_ccme

   extracted, blocked_core = spectra_ccme(topology_toy_model(1), ["r6"],
                                          tol=1e-4,
                                          consistency_type="topology")
   print(len(extracted.reactions), blocked_core)

.. code-block:: text

   6 []

Under the steady state the same call reports ``r6`` as blocked and returns an
empty model, since nothing in that network can carry flux at ``n=1``.

Both of these results are asserted in the test suite.

Next
----

* :doc:`consistency` — the consistency check in general
* :doc:`extraction` — building a context-specific model
