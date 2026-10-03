Choosing a formulation
======================

``problem_type`` decides what the extraction actually optimises. The choice
changes the model you get, and it changes what a weight means, so the two go
together.

.. list-table::
   :header-rows: 1
   :widths: 20 42 18 20

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
   * - ``"tradeOff"``
     - maximise weighted count of kept reactions
     - **any real**
     - 1 MILP
   * - ``"growthOptim"``
     - maximise objective flux minus weighted flux
     - non-negative
     - 1 LP

Seeing the difference
---------------------

.. code-block:: python

   from cobra.io import load_model
   from spectra_cobra import spectra_cc, spectra_me

   model = load_model("textbook")
   consistent = spectra_cc(model, tol=1e-3)

   for problem_type in ("minNetLP", "minNetMILP", "growthOptim", "tradeOff"):
       extracted = spectra_me(consistent, ["PGI"], tol=1e-3,
                              problem_type=problem_type, seed=0)
       print(f"{problem_type:<12} {len(extracted.reactions):>3} reactions")

.. code-block:: text

   minNetLP      27 reactions
   minNetMILP    16 reactions
   growthOptim   27 reactions
   tradeOff      87 reactions

Four very different models from the same core set. The spread is the point:
these are not four implementations of one idea.

minNetLP, the default
---------------------

Minimises :math:`\sum_i w_i t_i` subject to :math:`t_i \ge |v_i|` — an L1
relaxation of minimising the reaction count. One LP, so it is by far the
cheapest and the only one that is practical to run hundreds of times.

The catch is in the relaxation. L1 favours many small fluxes over a few large
ones, so it tends to keep more reactions than strictly necessary: 27 here
against the true minimum of 16.

minNetMILP, the exact one
-------------------------

Gives each candidate reaction a binary :math:`z_i` tied to its flux by
:math:`lb_i z_i \le v_i \le ub_i z_i`, so :math:`z_i = 0` pins the flux to
zero, and minimises :math:`\sum_i w_i z_i`. That is the reaction count itself,
not a relaxation, which is why it finds 16 where ``minNetLP`` finds 27.

It is a MILP. On the textbook model that is instant; on a genome-scale model
it does not finish at all, so use Gurobi or CPLEX and cap the effort with
``time_limit`` (seconds, default 7200). A solution found before the cap is
returned even if not proven optimal, with a warning. That is the normal
outcome at genome scale, and the incumbent is still a far smaller model than
``minNetLP`` gives — see :doc:`recon3d` for the measurements.

.. code-block:: python

   extracted = spectra_me(consistent, ["PGI"], tol=1e-3,
                          problem_type="minNetMILP", time_limit=600)

tradeOff, for evidence that cuts both ways
------------------------------------------

The others treat a weight as a cost, so it must be non-negative. ``tradeOff``
treats it as *evidence* and takes any real number: positive pushes a reaction
in, negative pushes it out. That is what you want when your omics data says
some reactions are likely present and others likely absent.

.. code-block:: python

   weights = {r.id: 1.0 for r in consistent.reactions}
   for rxn_id in ("SUCDi", "FRD7", "FUM", "MDH"):
       weights[rxn_id] = -5.0

   extracted = spectra_me(consistent, ["PGI"], tol=1e-3, weights=weights,
                          problem_type="tradeOff", seed=0)
   kept = {r.id for r in extracted.reactions}
   print(len(kept), [r for r in ("SUCDi", "FRD7", "FUM", "MDH") if r in kept])

.. code-block:: text

   83 []

All four penalised reactions are gone, and the rest stay because their +1
weights reward inclusion.

That last part explains the 87 in the table above: with the default weights of
+1 everywhere, ``tradeOff`` is asked to include as much as possible, so it
keeps essentially the whole model. **``tradeOff`` with default weights is not
a useful call** — it only does something meaningful once the weights carry
both signs.

growthOptim, for growth
-----------------------

Maximises :math:`\sum_{i \in c} w_i v_i - \sum_{i \notin c} w_i t_i`, where
:math:`c` are the reactions in ``model.objective``: reward the objective flux,
charge for everything else. Use it when the extracted model needs to grow,
not merely to be feasible.

The reward on an objective reaction is its *weight*, not its coefficient in
``model.objective`` — the objective is used only to find which reactions those
are, so the weight is what controls the trade-off. And unlike the other
formulations it charges for the absolute flux of *every* non-objective
reaction, core reactions included, not just the ones it is free to drop. A
model with no objective raises :class:`~spectra_cobra.SpectraError`.

.. warning::

   **The objective weight has to be scaled against the model.** The reward is
   ``weight x objective flux`` against a penalty of one per unit of flux
   everywhere else, and a genome-scale network at full growth carries
   *hundreds* of units of total flux. If the weight is too small, the
   optimum is simply not to grow at all:

   .. code-block:: python

      for weight in (1.0, 100.0, 1000.0, 5000.0):
          weights = {r.id: 1.0 for r in consistent.reactions}
          weights["Biomass_Ecoli_core"] = weight
          extracted = spectra_me(consistent, ["PGI"], tol=1e-3,
                                 weights=weights,
                                 problem_type="growthOptim", seed=0)
          print(f"{weight:>7g} {len(extracted.reactions):>3} reactions, "
                f"growth {extracted.slim_optimize():.4f}")

   .. code-block:: text

         1  27 reactions, growth 0.0000
       100  27 reactions, growth 0.0000
      1000  48 reactions, growth 0.8739
      5000  48 reactions, growth 0.8739

   Below roughly 1000 the biomass reaction is not even kept and the result is
   indistinguishable from ``minNetLP``; above it, growth jumps straight to the
   unconstrained maximum of 0.8739 and the model gains the 21 reactions needed
   to support it. **Always check that the extracted model actually grows**
   rather than assuming the formulation arranged it:

   .. code-block:: python

      assert extracted.slim_optimize() > 0, "growthOptim produced a dead model"

   MATLAB's example sets this weight to 10, which works there only because its
   toy model carries fluxes of order one.

Which one should I use?
-----------------------

* **Start with ``minNetLP``.** It is the default, it is one LP, and for most
  purposes a slightly larger model is not a problem.
* **Use ``minNetMILP``** when the model size is the result you care about, and
  you can afford a MILP.
* **Use ``tradeOff``** when you have signed evidence. Do not use it with
  default weights.
* **Use ``growthOptim``** when the extracted model has to grow.

A note on minNetDC
------------------

MATLAB offers a third route to the minimum-count objective, ``minNetDC``,
which delegates to the COBRA Toolbox's ``optimizeCardinality``. This package
does not provide it: cobrapy has no equivalent, and ``minNetMILP`` targets the
same objective and solves it exactly.

Next
----

* :doc:`alternative_solutions` — several models instead of one
* :doc:`recon3d` — what these cost at genome scale
