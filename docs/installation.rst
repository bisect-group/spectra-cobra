Installation
============

Requirements
------------

* Python 3.9 or newer
* `cobrapy <https://github.com/opencobra/cobrapy>`_ 0.29 or newer
* NumPy 1.20 or newer, and optlang 1.8 or newer (both come with cobrapy)

An LP solver is needed, and cobrapy ships with GLPK, so the consistency check
and the LP formulations work out of the box.

Download and installation
-------------------------

.. code-block:: bash

   git clone https://github.com/bisect-group/spectra-cobra.git
   cd spectra-cobra
   pip install -e .

Solvers
-------

Which solver you want depends on which formulations you plan to use.

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Function
     - Problem
     - Solver
   * - :func:`~spectra_cobra.spectra_cc`
     - LP
     - any, GLPK included
   * - ``minNetLP``, ``growthOptim``
     - LP
     - any, GLPK included
   * - ``minNetMILP``, ``tradeOff``
     - MILP
     - GLPK works; Gurobi or CPLEX much faster

The mixed-integer formulations solve on GLPK but can be slow on a
genome-scale model. To install a commercial solver through this package's
extras:

.. code-block:: bash

   pip install -e ".[gurobi]"   # or ".[cplex]"

Both offer free academic licences. Set the solver on the model before
calling anything:

.. code-block:: python

   model.solver = "gurobi"

Checking the installation
-------------------------

.. code-block:: bash

   pip install -e ".[development]"
   pytest

The suite runs every test against each solver cobrapy can see, so installing
Gurobi or CPLEX widens the coverage automatically.

.. code-block:: python

   >>> import spectra_cobra
   >>> spectra_cobra.__version__
   '0.1.0.dev0'
   >>> spectra_cobra.PROBLEM_TYPES
   ('minNetLP', 'minNetMILP', 'growthOptim', 'tradeOff')
