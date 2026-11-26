Formulations
============

The network inference problems :doc:`spectra_me` and :doc:`spectra_ccme` can
use, selected with ``problem_type``. Each is also callable on its own if you
want to supply the reaction directions yourself.

.. list-table::
   :header-rows: 1
   :widths: 22 18 40 20

   * - ``problem_type``
     - Function
     - Objective
     - Weights
   * - ``"minNetLP"``
     - :func:`~spectra_cobra.min_net_lp`
     - minimise weighted total absolute flux
     - non-negative
   * - ``"minNetMILP"``
     - :func:`~spectra_cobra.min_net_milp`
     - minimise weighted reaction count
     - non-negative
   * - ``"tradeOff"``
     - :func:`~spectra_cobra.trade_off`
     - maximise net weight: positive in, negative out
     - **any real**
   * - ``"optimBiomass"``
     - :func:`~spectra_cobra.optim_biomass`
     - maximise objective flux minus weighted flux
     - non-negative

See :doc:`../tutorials/formulations` for which to choose and what the weight
conventions mean in practice.

Indicators, and sharing one
---------------------------

The mixed-integer formulations decide membership with a binary per free
reaction. Two arguments change that.

``indicator_reactions``
   Restrict the binaries to these reactions. Everything else free stays
   continuous and is kept whatever it does. One binary per *organism* rather
   than per reaction is what makes
   :func:`~spectra_cobra.minimal_microbiome` tractable on a community.

``indicator_groups``
   Give the reactions of each group **one shared binary**, so they are kept
   or dropped together and cost (or score) the group's weight once however
   many of them there are. Supported by both ``minNetMILP`` and ``tradeOff``,
   and mutually exclusive with ``indicator_reactions``.

.. math::

   \text{minimize} \quad \sum_g w_g z_g
   \qquad \text{s.t.} \quad
   lb_i\, z_g \le v_i \le ub_i\, z_g \quad \text{for every } i \in g

Sharing a binary is how a requirement that holds in several conditions at
once is expressed: replicate the network once per condition, group each
reaction's copies together, and the objective counts reactions while the
constraints are satisfied separately in every condition. That is exactly
what :func:`~spectra_cobra.minimal_reactome` does.

With neither argument each group has a single member keyed by the reaction's
own identifier, which is the ordinary formulation and why
:attr:`~spectra_cobra.MilpSolution.selected` reads as a set of reactions.
With groups it is a set of group keys.

In ``tradeOff``, switching a group binary :math:`z_g = 1` requires every
governed member :math:`i \in g` to carry at least ``tol`` of flux (with
:math:`a_i + b_i = z_g` on reversible members).

Directions
----------

Called directly, each function takes a ``directions`` dict mapping every
reaction identifier to one of:

``1``
   The reaction must carry positive oriented flux, of at least ``tol``.
``-1``
   It must carry negative oriented flux, of at least ``tol`` in magnitude.
``0``
   It is free; the formulation decides whether to keep it.

The reactions with a nonzero direction are the core set.
:func:`~spectra_cobra.spectra_me` builds this dict for you from the flux its
first phase accumulates, so you only need it if you are driving the
formulations yourself.

.. code-block:: python

   from spectra_cobra import min_net_lp

   directions = {r.id: 0 for r in model.reactions}
   directions["PGI"] = 1          # force positive flux through PGI
   weights = {r.id: 1.0 for r in model.reactions}

   kept = min_net_lp(model, directions, weights, tol=1e-3)

Note that these return a **set of reaction identifiers**, not a model —
:func:`~spectra_cobra.spectra_me` is what turns that into one.

The mathematics
---------------

**minNetLP** computes "With-Chain" median-normalised non-negative weights
(scaling positive open non-core reaction weights by their degree-2 linear
chain length in ``stoichiometry`` mode) and solves a warm-started
Split-Variable iteratively reweighted L1 sequence
(``deltas = [1e-3, 1e-5, 1e-6, 1e-8]``) followed by tolerance-derived support
certification:

.. math::

   \text{minimize} \quad \sum_i w_i^{(t)} |v_i|
   \qquad \text{s.t.} \quad S v = 0, \; v \in B

**minNetMILP** minimises the count itself, with a binary per free reaction:

.. math::

   \text{minimize} \quad \sum_i w_i z_i
   \qquad \text{s.t.} \quad lb_i z_i \le v_i \le ub_i z_i, \;
                           z_i \in \{0, 1\}

**tradeOff** maximises the net weight of the included reactions, so a
positive weight pushes a reaction in and a negative one pushes it out. A free
irreversible reaction gets one binary; a free *reversible* one gets two more,
:math:`a_i` for running forward and :math:`b_i` for backward, with
:math:`a_i + b_i = z_i` so at most one direction is active:

.. math::

   \text{maximize} \quad & \sum_i w_i z_i \\
   \text{s.t.} \quad & \hat v_i \ge \varepsilon a_i + lb_i b_i \\
                     & \hat v_i \le ub_i a_i - \varepsilon b_i

where :math:`\hat v` is the oriented flux. With
:math:`a_i = b_i = 0` these pin the flux to zero. Note that because
``tradeOff`` maximises :math:`\sum_i w_i z_i` subject only to steady-state
mass balance and bounds, any disconnected internal stoichiometrically balanced
cycle (Type III loop) with net positive weight can circulate flux of at least
:math:`\varepsilon` independently of the core set and boundary exchanges. In
both ``minNetMILP`` and ``tradeOff``, the returned reaction set is the union of
core reactions, switched-on indicators (:math:`z_i \ge 0.5`), and reactions
with :math:`|v_i| \ge 10^{-7}\varepsilon`, while ``previous_solutions``
exclusion cuts constrain the switched-on indicator sets directly.

**optimBiomass** rewards the objective reactions :math:`c` and charges for
everything else:

.. math::

   \text{maximize} \quad \sum_{i \in c} w_i v_i - \sum_{i \notin c} w_i t_i
   \qquad \text{s.t.} \quad t_i \ge |v_i|

Reference
---------

.. autofunction:: spectra_cobra.min_net_lp

.. autofunction:: spectra_cobra.min_net_milp

.. autofunction:: spectra_cobra.trade_off

.. autofunction:: spectra_cobra.optim_biomass

See also
--------

* :doc:`../tutorials/formulations` — choosing between them, with examples
* :doc:`spectra_me` — the usual way to invoke them
