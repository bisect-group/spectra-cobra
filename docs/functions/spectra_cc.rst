spectra_cc
==========

Removes the blocked reactions of a metabolic network, leaving a flux
consistent model.

Usage
-----

.. code-block:: python

   from spectra_cobra import spectra_cc

   consistent = spectra_cc(model, tol=1e-4)

How it works
------------

Two LPs alternate, each driving as many of the still-unexplained reactions as
possible towards the flux threshold at once, until that set stops shrinking:

.. math::

   \text{maximize}   &\quad \sum_{i \in J} w_i z_i \\
   \text{subject to} &\quad z_i \le \varepsilon, \;
                            v_i \ge z_i \quad \forall i \in J \\
                     &\quad S v = 0, \; v \in B

.. math::

   \text{minimize}   &\quad \sum_{i \in J^{rev}} w_i z_i \\
   \text{subject to} &\quad z_i \ge -\varepsilon, \;
                            v_i \le z_i \quad \forall i \in J^{rev} \\
                     &\quad S v = 0, \; v \in B

:math:`J` is the set of reactions not yet shown to carry flux,
:math:`J^{rev}` the reversible ones among them, and :math:`w_i` are drawn
uniformly from [1, 1.1] to break ties.

No reaction flipping is needed. The forward LP's auxiliary variable
:math:`z_i` is unbounded from below, so a reaction that can only carry
negative flux simply fails to contribute to that LP and is picked up by the
reverse LP instead.

Reference
---------

.. autofunction:: spectra_cobra.spectra_cc

.. autofunction:: spectra_cobra.consistent_reaction_ids

.. autofunction:: spectra_cobra.blocked_reaction_ids

See also
--------

* :doc:`../tutorials/consistency` — worked examples
* :doc:`../tutorials/topology` — the two consistency types compared
* :doc:`spectra_ccme` — consistency checking folded into extraction
