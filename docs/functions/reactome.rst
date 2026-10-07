minimal_reactome
================

Cuts a model down to the smallest set of reactions that still meets a stated
requirement: a medium, a growth rate, the production of a metabolite, a list
of metabolic tasks, or any combination.

:doc:`../tutorials/minimal_reactome` works all of this through, from the
textbook *E. coli* core model to iJO1366.

Usage
-----

.. code-block:: python

   from cobra.io import load_model
   from spectra_cobra import minimal_reactome

   model = load_model("textbook")
   model.tolerance = 1e-9

   result = minimal_reactome(
       model,
       medium={"EX_glc__D_e": -10.0, "EX_o2_e": -1000.0,
               "EX_nh4_e": -1000.0, "EX_pi_e": -1000.0,
               "EX_h2o_e": -1000.0, "EX_h_e": -1000.0,
               "EX_co2_e": -1000.0},
       growth_fraction=0.5,
   )
   print(result.summary())
   smaller = result.model

The model handed in is not modified. The result is verified against every
requirement before it is returned, and ``result.unsatisfied`` names anything
it misses.

How it works
------------

Minimising a count is a mixed-integer problem. Every reaction gets a binary
:math:`y_j` tied to its flux, so that switching the binary off pins the flux
to zero:

.. math::

   \text{minimize} \quad \sum_j w_j y_j
   \qquad \text{s.t.} \quad
   S v = 0, \quad lb_j\, y_j \le v_j \le ub_j\, y_j

This is :func:`~spectra_cobra.min_net_milp` with no core reactions: the
requirements enter as bounds rather than as a core set, which is what lets
them be a growth rate or a secretion rate rather than merely "must carry some
flux".

A metabolic task defines its own medium, so it cannot share a flux vector
with the growth requirement. Tasks therefore become separate **conditions**,
and so does each medium when more than one is given: the network is
replicated once per condition, and one binary per reaction governs every
copy at once.

.. math::

   \text{minimize} \quad \sum_j w_j y_j
   \qquad \text{s.t.} \quad
   S^{(k)} v^{(k)} = 0, \quad
   lb^{(k)}_j\, y_j \le v^{(k)}_j \le ub^{(k)}_j\, y_j
   \quad \forall k

Solving the conditions jointly is what makes the answer minimal across them.
On the textbook model with one growth condition and three tasks, the joint
solve keeps 51 reactions where the union of four separate solves keeps 58.

Keeping reactions
-----------------

``keep_reactions`` is a **core set**, not merely a retain list: each reaction
is given a direction by the same two LPs :func:`~spectra_cobra.spectra_me`
uses, and then forced to carry at least ``tol``. This applies in the growth
conditions only — a task closes every exchange, so forcing a boundary
reaction inside one would be infeasible by construction. Anything that cannot
reach ``tol`` in a growth condition is retained anyway and reported in
``result.kept_but_blocked``.

Preprocessing
-------------

``preprocess`` is on by default and settles two things exactly before the
MILP starts. Reactions that can carry no flux in **any** condition are
deleted (:func:`~spectra_cobra.spectra_cc` per condition, intersected);
reactions **some** condition is infeasible without are forced in and given
no binary (unioned). On iJO1366 on glucose that is 2057 blocked and 407
essential, leaving 119 binaries out of 2583 — and the objective is 28
either way, so the optimum is untouched.

It costs LPs, so on a problem that already closes quickly it is a net loss:
7 s against 14 s for that same glucose run. Turn it off there.

``result.essential`` reports the reactions proved to be in every possible
answer.

Related
-------

* :func:`~spectra_cobra.gapfill_for_tasks` — the opposite operation: adds
  reactions until the tasks pass, rather than removing them while they do
* :func:`~spectra_cobra.minimal_microbiome` — the same question asked of
  organisms instead of reactions
* :func:`~spectra_cobra.spectra_me` — extraction around a core set, when you
  know which reactions must be present rather than what the network must do

Reference
---------

.. autofunction:: spectra_cobra.minimal_reactome

.. autoclass:: spectra_cobra.MinimalReactome
   :members:

.. autodata:: spectra_cobra.reactome.DEFAULT_GROWTH_FRACTION

.. autodata:: spectra_cobra.reactome.INCLUSION_CUTOFF

See also
--------

* :doc:`../tutorials/minimal_reactome` — the worked walkthrough
* :doc:`../tutorials/formulations` — what ``minNetMILP`` is doing
* :doc:`tasks` — stating what a network must be able to do
