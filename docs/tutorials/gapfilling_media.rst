Gap-filling towards growth in a medium
======================================

:doc:`gapfilling_tasks` fills gaps so a model can perform a list of tasks.
This page does the same thing with a different requirement: the model must
*grow*, in each of a panel of media.

That is the setting of a Biolog-style phenotyping experiment. You know which
carbon sources the organism grows on, your draft reconstruction fails on some
of them, and you want the smallest set of reactions from a universal model
that accounts for the difference.

The shape of the problem
------------------------

For each medium:

1. constrain the **universal** model to that medium,
2. find the reactions essential for growth in it,
3. ask for the cheapest set of additions that lets the draft grow,

and union the additions over the media. The essential reactions are made free
rather than forced: they cost nothing to include, so the solver never routes
around one only to pay for an alternative.

The core set is just the biomass reaction. Everything else follows from
requiring it to carry flux.

.. important::

   **Demand growth with a bound, not with tol.** ``min_growth`` puts a lower
   bound on the biomass reaction; ``tol`` stays at its usual extraction
   value. It is tempting to make biomass the core reaction and set ``tol``
   to the growth rate instead, and it appears to work — but ``tol`` is the
   flux *every* core reaction must carry, so that silently demands the same
   rate of any other core reaction you supply. The last section shows the
   difference.

A worked example
----------------

On the textbook model, across nine carbon sources:

.. code-block:: python

   import numpy as np
   from cobra.flux_analysis import single_reaction_deletion
   from cobra.io import load_model
   from spectra_cobra import spectra_me

   universal = load_model("textbook")
   BIOMASS = "Biomass_Ecoli_core"
   SALTS = ["EX_pi_e", "EX_h2o_e", "EX_h_e", "EX_nh4_e", "EX_o2_e", "EX_co2_e"]
   CARBON = ["EX_glc__D_e", "EX_ac_e", "EX_succ_e", "EX_mal__L_e", "EX_fum_e",
             "EX_pyr_e", "EX_akg_e", "EX_gln__L_e", "EX_glu__L_e"]

   def medium(model, carbon, uptake=10.0):
       """Close every exchange, then open the salts and one carbon source."""
       for reaction in model.boundary:
           reaction.lower_bound = 0.0
       for ex_id in SALTS:
           model.reactions.get_by_id(ex_id).lower_bound = -1000.0
       model.reactions.get_by_id(carbon).lower_bound = -uptake

Closing everything first matters: a medium is defined by what it *lacks* as
much as by what it supplies, and a model's default bounds usually leave far
more open than you intend.

**Reactions essential for growth**, by single deletion in each medium:

.. code-block:: python

   essential = set()
   for carbon in usable:
       with universal:
           medium(universal, carbon)
           result = single_reaction_deletion(universal, processes=1)
       ids = [next(iter(s)) for s in result["ids"]]
       essential |= {i for i, g in zip(ids, result["growth"].to_numpy())
                     if np.isnan(g) or g < 0.01}

.. code-block:: text

   EX_glc__D_e      18 essential
   EX_ac_e          37 essential
   EX_succ_e        30 essential
   ...
   union: 54

Acetate needs twice as many reactions as glucose, because growing on a
two-carbon compound requires the glyoxylate shunt and gluconeogenesis that
glucose makes unnecessary. Each medium contributes its own constraints, which
is the reason to gap-fill against a panel rather than a single condition.

**A draft with gaps**, and it grows nowhere:

.. code-block:: text

   removed ['ENO', 'ACONTb', 'SUCDi', 'NH4t']
   media the draft cannot grow in: 9/9

**Gap-fill**, medium by medium, letting the draft, the essential reactions,
the exchanges and anything already added be free:

.. code-block:: python

   added = set()
   for carbon in usable:
       with universal:
           medium(universal, carbon)
           free = (draft_ids | essential | added
                   | {rxn.id for rxn in universal.boundary})
           weights = {rxn.id: (0.0 if rxn.id in free else 1.0)
                      for rxn in universal.reactions}
           filled = spectra_me(universal, [BIOMASS], tol=0.1, weights=weights,
                               problem_type="minNetMILP", seed=0)
       added |= {rxn.id for rxn in filled.reactions} - draft_ids - added

   # The gap-filled model is the draft plus the additions.
   gapfilled = universal.copy()
   gapfilled.remove_reactions(
       sorted({rxn.id for rxn in universal.reactions} - draft_ids - added),
       remove_orphans=False)

.. code-block:: text

   EX_glc__D_e      added ['ACONTb', 'ENO', 'NH4t', 'SUCDi']

   added in total: ['ACONTb', 'ENO', 'NH4t', 'SUCDi']
   grows in 9/9 media

All four holes were found from the glucose condition alone, and the remaining
eight media needed nothing further. That is the usual pattern: the first
medium does most of the work and the rest confirm it, which is why
``added`` is carried into the next iteration's free set.

.. note::

   Adding to ``free`` as you go is what keeps the result small. Without it
   each medium prices the previous medium's additions all over again, and
   the union grows with reactions that duplicate work already paid for.

Why the essential reactions are free, not core
----------------------------------------------

In :doc:`gapfilling_tasks` the essential reactions become the *core set*:
they are forced to carry flux. Here they are given **weight zero** instead,
which only makes them free to include.

The difference is what is being required. There, the task is the requirement
and the essential reactions are how it is expressed. Here the requirement is
growth, stated directly by making biomass the core, and the essential
reactions are a hint to the solver: these are reactions growth is known to
need in some condition, so do not charge for them. Forcing all of them to
carry flux at once would be wrong — they come from different media, and no
single flux distribution runs all of them.

At genome scale
---------------

The same workflow on iJO1366 (2583 reactions) across twelve M9 carbon
sources — glucose, acetate, succinate, glycerol, galactose, malate,
fumarate, pyruvate, ribose, xylose, fructose and lactate — with a draft made
by deleting 150 random non-exchange reactions:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * -
     - Reactions
     - Added
     - Media grown in
   * - universal
     - 2583
     - —
     - 12 / 12
   * - draft
     - 2433
     - —
     - **0 / 12**
   * - gap-filled (``minNetMILP``)
     - 2435
     - **2**
     - **12 / 12**

Two reactions. The draft was missing 150 and could not grow on anything, and
two well-chosen additions restored every condition.

That is worth sitting with, because it is the central caveat of gap-filling:
the result is a *minimal repair*, not a reconstruction of what was removed.
148 of the deleted reactions were never missed, because the network routed
around them. Gap-filling tells you what the model needs in order to grow,
which is not the same question as what the organism actually does.

Finding the essential reactions took 45 s (about 290 per medium, 331 in
union) and the twelve gap-fills 54 s in total.

Why not tradeOff here
---------------------

``tradeOff`` is the natural-looking choice, since signed weights express
"keep the draft, charge for additions" directly. It does not work for this
problem, and the reason is structural rather than a matter of tuning: on all
twelve media it came back **infeasible**.

``tradeOff`` ties each candidate reaction to a binary with

.. math::

   \varepsilon z_i \le \hat v_i \le ub_i z_i

so a reaction that is included must carry at least ``tol``. Growth does not
work like that. In a full-growth solution on iJO1366, 442 reactions carry
flux, but 232 of them carry less than 0.1, 24 carry less than ``1e-4``, and
the smallest carries ``1.96e-06`` — trace demands for cofactors and
minerals that are nonetheless structurally necessary. Requiring every one of
them to reach ``tol`` has no solution.

The prediction that follows is exact: ``tradeOff`` should become feasible
only once ``tol`` drops below that smallest flux, and it does. At
``tol=1e-5`` it is infeasible; at ``tol=1e-6`` it returns a model. By then
``tol`` is far below any growth rate worth asking for, so the formulation has
not been rescued, only emptied of meaning.

``minNetMILP`` has no such requirement. Its constraint is
:math:`lb_i z_i \le v_i \le ub_i z_i`, which pins the flux to zero when the
reaction is out but never imposes a floor when it is in, so a trace flux is
perfectly acceptable.

Use ``tradeOff`` where the requirement is a set of reactions rather than
biomass — :doc:`gapfilling_tasks` is that case, and it works there — and
``minNetMILP`` whenever growth is what you are asking for.

Growth is a bound, not a tolerance
----------------------------------

:func:`~spectra_cobra.gapfill_for_growth` demands growth through
``min_growth``, a lower bound on the biomass reaction, and leaves ``tol`` at
the value any extraction would use. The alternative — making biomass the
core reaction and raising ``tol`` to the growth rate — looks equivalent and
is not, for two reasons.

**It spreads to every core reaction.** ``tol`` is the flux each core
reaction must carry. Raise it to 0.1 to mean "grow at 0.1" and you have also
demanded 0.1 through every other reaction you named as core, which is rarely
what you meant and often infeasible.

**It makes the answer move with the tolerance.** With the bound in place the
result stops depending on ``tol`` at all. On iJO1366, across a glucose
minimal medium:

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - ``tol``
     - Reactions added
     - Growth
     - ``tol`` as the demand
   * - ``1e-4``
     - 15
     - 0.683
     - **0.000** (silently broken)
   * - ``1e-5``
     - 15
     - 0.683
     - 0.623
   * - ``1e-6``
     - 15
     - 0.683
     - —
   * - ``1e-7``
     - 15
     - 0.683
     - —

The left column is the bound doing the work: identical at every tolerance.
The right column is what the same problem gives when ``tol`` carries the
demand instead, and the ``1e-4`` row is the trap — the LP pins growth at
``tol``, every supporting flux shrinks with it, the smallest fall below the
cutoff at which a reaction counts as used, and the chains they held up
break. Nothing errors, and the broken model is the *smaller* one.

.. note::

   ``min_growth`` defaults to 0.1, which is a plausible rate rather than a
   measured one. If you have an experimental growth rate for the organism
   and condition, pass it: the point of the default is to be a reasonable
   demand, not to be right about your organism.

Next
----

* :doc:`gapfilling_tasks` — the same idea with a task list as the requirement
* :doc:`formulations` — what the objectives optimise
* :doc:`recon3d` — tolerances and cutoffs at genome scale
