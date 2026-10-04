Differences from the MATLAB implementation
==========================================

The optimisation problems are the same in both implementations. The machinery
around them is built on cobrapy's solver interface rather than on raw
matrices, and a few behaviours are deliberately different.

Nothing mutates your model
--------------------------

MATLAB begins every routine by negating the stoichiometric columns of
reactions that can only carry negative flux, and swapping their bounds:

.. code-block:: matlab

   IrR = model.ub <= 0;
   model.S(:, IrR) = -model.S(:, IrR);
   model.ub(IrR) = -model.lb(IrR);
   model.lb(IrR) = -temp;

That rewrite gives the algorithm a canonical "try this way first" direction.
Here the same thing is a per-reaction sign applied to the reaction's flux
expression, so the input model is never touched. The two are equivalent
because negating a column and negating its variable cancel:

.. math::

   \sum_i (\pm S_{:,i})(\pm v_i) = \sum_i S_{:,i} v_i

which also shows why the mass balance needs no sign handling, and only the
terms that single out a direction do.

Failures raise instead of returning NaN
---------------------------------------

MATLAB warns and returns ``NaN`` or an empty solution when an LP does not
reach optimality. Its callers then misread that. In ``spectraME``,
``abs(NaN) >= tol`` is false, so the loop waiting for core reactions to be
explained never shrinks its working set and **never terminates**.

Here a non-optimal solve raises, so an impossible core set gives you
:class:`~spectra_cobra.SpectraInfeasibleCoreError` with an explanation rather
than a hang.

Mixed-integer results are read off the binaries
-----------------------------------------------

MATLAB derives the extracted reaction set from the flux,
``abs(x) >= tol * 1e-7``, for the MILP formulations as well as the LP ones. At
a default ``tol`` that cutoff is around ``1e-10``, which is far below any
solver's integrality tolerance: a binary resting at ``1e-7`` against a flux
bound of 1000 leaks ``1e-4`` of apparent flux and reads as a kept reaction.

That matters most for pathway exclusion, whose constraint is written over
those same binaries. If the recorded reaction set disagrees with them, the
constraint can be satisfied without changing anything, and **the same model
comes back while being reported as an alternative**. This port reads the
mixed-integer result off the binaries directly, and floors the LP formulations'
flux cutoff at ``model.tolerance``.

MATLAB partly masks the problem with
``changeCobraSolverParams('MILP', 'feasTol', 1e-9)``, which does not carry
over to optlang.

Fixed MATLAB bugs
-----------------

``spectraCCME.m`` passes a misspelled ``steadyState`` to its ``growthOptim``
and ``tradeOff`` branches where the variable is ``steadystate``. MATLAB is
case-sensitive, so in MATLAB those two problem types raise an
undefined-variable error rather than running. Both work here.

Not ported
----------

``minNetDC``
   MATLAB delegates it to the COBRA Toolbox's ``optimizeCardinality``, which
   cobrapy has no equivalent of. ``minNetMILP`` targets the same objective and
   solves it exactly.

``spectraME2.m``
   An older, LP-only subset of ``spectraME``, which covers everything it does.

Smaller changes
---------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Area
     - Change
   * - Reproducibility
     - Every routine takes ``seed``, so the randomised objective
       coefficients can be pinned, in a way that holds between processes as
       well as within one. MATLAB's ``unifrnd`` cannot be pinned at all.
   * - Pathway exclusion
     - Returns the solutions that exist when the network runs out, like
       MATLAB's ``break``, rather than raising.
   * - Validation
     - An unknown ``consistency_type``, ``problem_type``, core reaction or
       weighted reaction raises immediately. MATLAB leaves ``steadystate``
       undefined on a bad ``consType``, failing later and less clearly.
   * - Weights and cores
     - Given as dicts and iterables keyed by reaction identifier, not as
       index-aligned vectors, so they cannot silently misalign with the
       model.

Verified parity
---------------

``tests/test_matlab_parity.py`` rebuilds the MATLAB toy-model experiments and
asserts this port reaches the same answers, for
``SPECTRA_CC_topology_vs_stoichiometry.m``,
``SPECTRA_ME__topology_vs_stoichiometry.m`` and
``Objective_diff_toy_models.m``. The expected values are derived from each
model's own stoichiometry rather than copied from a MATLAB run, so they are
checked rather than assumed. See :doc:`tutorials/topology` for the most
instructive of them.
