flux_reducer
============

Finds a sparse flux distribution for a model, by minimising the sum of
absolute fluxes.

Usage
-----

.. code-block:: python

   from spectra_cobra import flux_reducer

   fluxes = flux_reducer(model)
   active = {r: v for r, v in fluxes.items() if abs(v) > 1e-9}

Returns a dict mapping reaction identifier to flux.

How it works
------------

.. math::

   \text{minimize} \quad \sum_{i \notin rev} \hat v_i + \sum_{i \in rev} t_i
   \qquad \text{s.t.} \quad t_i \ge |v_i|, \; S v = 0, \; v \in B

An irreversible reaction is charged for its own oriented flux
:math:`\hat v_i`, which is non-negative and so already an absolute value.
Only a reversible reaction needs an auxiliary variable, which spares the
problem one variable per irreversible reaction.

Reference
---------

.. autofunction:: spectra_cobra.flux_reducer
