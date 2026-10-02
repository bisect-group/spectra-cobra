Checking consistency
====================

A reaction is *blocked* if it cannot carry flux in any steady state. A
universal model assembled from a database usually has plenty of them, and they
matter: :func:`~spectra_cobra.spectra_me` requires a consistent model, and a
blocked reaction in your core set makes the extraction infeasible.

Removing the blocked reactions
------------------------------

.. code-block:: python

   from cobra.io import load_model
   from spectra_cobra import spectra_cc

   model = load_model("textbook")
   print(len(model.reactions))

   consistent = spectra_cc(model, tol=1e-3)
   print(len(consistent.reactions))

.. code-block:: text

   95
   87

Eight of the textbook model's 95 reactions cannot carry flux. The returned
model is a copy with those removed, and with any metabolite left without a
reaction removed too; the input model is untouched.

Which reactions were blocked?
-----------------------------

:func:`~spectra_cobra.blocked_reaction_ids` answers that directly, without
building a model:

.. code-block:: python

   from spectra_cobra import blocked_reaction_ids

   print(blocked_reaction_ids(model, tol=1e-3))

.. code-block:: text

   ['EX_fru_e', 'EX_fum_e', 'EX_gln__L_e', 'EX_mal__L_e',
    'FRUpts2', 'FUMt2_2', 'GLNabc', 'MALt2_2']

These are transporters and exchanges for metabolites the textbook model's
default medium does not supply, paired with the uptake reactions that depend
on them.

Counting the work
-----------------

:func:`~spectra_cobra.consistent_reaction_ids` returns the identifiers and the
number of LPs that were solved:

.. code-block:: python

   from spectra_cobra import consistent_reaction_ids

   ids, n_lps = consistent_reaction_ids(model, tol=1e-3)
   print(len(ids), n_lps)

.. code-block:: text

   87 8

Four rounds of two LPs each. The count is a useful sanity check on a large
model: it grows very slowly with model size, so a run that takes many rounds
is worth a second look.

Agreement with cobrapy
----------------------

With the detection cutoff set to the solver tolerance, this agrees exactly
with cobrapy's own ``fastcc``:

.. code-block:: python

   from cobra.flux_analysis import fastcc

   spectra_ids, _ = consistent_reaction_ids(model, tol=1e-3,
                                            detection_cutoff=1e-7)
   fastcc_model = fastcc(model, 1e-3, 1e-7)

   print(spectra_ids == {r.id for r in fastcc_model.reactions})

.. code-block:: text

   True

The default ``detection_cutoff`` is ``0.99 * tol``, which is MATLAB's
criterion and a stricter one: it also discards reactions that can carry *some*
flux but never as much as ``tol``. Which you want depends on whether ``tol``
is a meaningful flux magnitude for your problem or just a numerical floor.

Choosing ``tol``
----------------

``tol`` is the flux threshold: the magnitude the LPs drive reactions towards,
and by default also the magnitude a reaction must reach to count as
consistent. Two things to keep in mind.

* It must sit comfortably above the solver's own tolerance, or you are
  measuring numerical noise. ``model.tolerance`` is ``1e-7`` by default, so
  ``tol`` of ``1e-4`` or larger is safe.
* Raising it tightens the definition. A reaction whose maximum attainable flux
  is ``1e-5`` is consistent at ``tol=1e-6`` and blocked at ``tol=1e-4``.

Reproducibility
---------------

The LPs use randomised objective coefficients to break ties between
reactions, so the flux distributions differ between runs even though the
consistent set does not. Pass ``seed`` to pin them:

.. code-block:: python

   ids_a, _ = consistent_reaction_ids(model, tol=1e-3, seed=0)
   ids_b, _ = consistent_reaction_ids(model, tol=1e-3, seed=0)
   assert ids_a == ids_b

Next
----

* :doc:`topology` — the accumulation condition, for topology-based gap filling
* :doc:`extraction` — building a context-specific model
