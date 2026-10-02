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

Both of these results are asserted in the test suite, in
``tests/test_matlab_parity.py``, against the MATLAB experiments they come
from.

Next
----

* :doc:`consistency` — the consistency check in general
* :doc:`extraction` — building a context-specific model
