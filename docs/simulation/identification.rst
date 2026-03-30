Parameter Identification
========================

Simcoon provides a Python-based parameter identification workflow using
``scipy.optimize`` and a generic key-based file templating system. This
replaces the former built-in C++ genetic algorithm and allows users to
choose any optimizer, cost function, or external simulation tool.

Overview
--------

The identification workflow consists of three components:

1. **Key system** (``Parameter`` / ``Constant``): generic file templating
   that replaces placeholders with parameter values in any input file
2. **Forward model**: simcoon solver, ``L_eff`` homogenization, or any
   external tool (e.g., fedoo for FEMU)
3. **Optimizer**: ``scipy.optimize.differential_evolution``, or any other
   scipy/scikit-learn optimizer

.. code-block:: none

   Template files (keys/)          Working files (data/)
   ┌──────────────────┐    copy    ┌──────────────────┐
   │  ... @2p ...     │ ────────>  │  ... @2p ...     │
   │  ... @3p ...     │            │  ... @3p ...     │
   └──────────────────┘            └──────────────────┘
                                          │ apply
                                          v
                                   ┌──────────────────┐
                                   │  ... 73000 ...   │  ──> Forward model
                                   │  ... 0.22 ...    │      (solver, L_eff, ...)
                                   └──────────────────┘

Key System
----------

The key system decouples the optimizer from the simulation tool. Template
files contain alphanumeric placeholders (keys) that are replaced at each
iteration with the current parameter values.

Parameter class
^^^^^^^^^^^^^^^

.. code-block:: python

   from simcoon.parameter import Parameter, read_parameters, copy_parameters, apply_parameters

   # Create parameters programmatically
   params = [
       Parameter(number=0, bounds=(100, 300), key="@E",
                 sim_input_files=["material.dat"]),
       Parameter(number=1, bounds=(0.1, 0.4), key="@nu",
                 sim_input_files=["material.dat"]),
   ]

   # Or read from a file
   params = read_parameters("data/parameters.inp")

**Parameter attributes:**

- ``number``: parameter index
- ``bounds``: ``(min, max)`` tuple — used as optimizer bounds
- ``key``: placeholder string in template files (e.g., ``@E``, ``@0p``)
- ``sim_input_files``: list of files containing this key
- ``value``: current value (defaults to midpoint of bounds)

**Parameters file format** (``parameters.inp``):

.. code-block:: none

   #Number  #min     #max     #key  #number_of_files  #files
   0        100      300      @E    1                  material.dat
   1        0.1      0.4      @nu   1                  material.dat

Constant class
^^^^^^^^^^^^^^

Constants are fixed values (not optimized) that also use the key system:

.. code-block:: python

   from simcoon.constant import Constant, read_constants, copy_constants, apply_constants

The ``Constant`` class is a ``NamedTuple`` with fields: ``number``, ``key``,
``input_values``, ``value``, ``sim_input_files``.

File operations
^^^^^^^^^^^^^^^

.. code-block:: python

   # 1. Copy template files from keys/ to data/
   copy_parameters(params, src_path="keys", dst_path="data")

   # 2. Replace keys with current values
   params[0].value = 200.0  # set by optimizer
   params[1].value = 0.3
   apply_parameters(params, dst_path="data")

This works with **any** file format — simcoon input files, FE meshes, JSON
configs, etc. The key system is deliberately simple: it performs string
replacement, making it compatible with any simulation tool.


Identification Workflow
-----------------------

Basic example
^^^^^^^^^^^^^

.. code-block:: python

   import numpy as np
   from scipy.optimize import differential_evolution
   import simcoon as sim
   from simcoon.parameter import Parameter, copy_parameters, apply_parameters

   # Known: experimental data
   c_exp = np.array([0.0, 0.1, 0.2, 0.3])
   E_exp = np.array([2250, 2580, 3390, 4480])

   # Define parameters with keys
   params = [
       Parameter(0, bounds=(0, 1), key="@c_m",
                 sim_input_files=["Nellipsoids0.dat"]),
       Parameter(1, bounds=(0, 1), key="@c_f",
                 sim_input_files=["Nellipsoids0.dat"]),
       Parameter(2, bounds=(10000, 200000), key="@Ef",
                 sim_input_files=["Nellipsoids0.dat"]),
   ]

   props = np.array([2, 0, 50, 50, 0], dtype="float")

   def cost_function(x):
       E_f = x[0]
       E_pred = np.zeros(len(c_exp))
       for i, c in enumerate(c_exp):
           params[0].value = 1 - c
           params[1].value = c
           params[2].value = E_f
           copy_parameters(params, "keys", "data")
           apply_parameters(params, "data")
           L = sim.L_eff("MIMTN", props, 0, 0., 0., 0.)
           E_pred[i] = sim.L_iso_props(L).flatten()[0]
       return np.mean((E_pred - E_exp) ** 2)

   result = differential_evolution(cost_function, bounds=[(10000, 200000)])
   print(f"Identified E_f = {result.x[0]:.0f} MPa")

Choosing an optimizer
^^^^^^^^^^^^^^^^^^^^^

``scipy.optimize`` provides several global and local optimizers:

- ``differential_evolution``: robust global optimizer, derivative-free,
  handles bounds naturally. Good default choice.
- ``dual_annealing``: simulated annealing, good for highly multimodal problems
- ``minimize`` with ``method='L-BFGS-B'``: fast local optimizer with bounds,
  use when you have a good initial guess

The ``polish=True`` option in ``differential_evolution`` automatically
refines the result with L-BFGS-B after the global search.

Cost function
^^^^^^^^^^^^^

The cost function is entirely user-defined. Common choices:

.. code-block:: python

   # Simple MSE
   mse = np.mean((y_pred - y_exp) ** 2)

   # Normalized MSE (for multi-experiment fitting)
   nmse = np.mean((y_pred - y_exp) ** 2) / np.var(y_exp)

   # Using sklearn for convenience
   from sklearn.metrics import mean_squared_error
   mse = mean_squared_error(y_exp, y_pred)


Using with External Solvers
----------------------------

The key system works with any simulation tool. For example, with fedoo
(Finite Element Model Updating — FEMU):

.. code-block:: python

   import subprocess
   from simcoon.parameter import Parameter, copy_parameters, apply_parameters

   params = [
       Parameter(0, bounds=(100e3, 300e3), key="@E",
                 sim_input_files=["material.json"]),
   ]

   def run_fedoo_and_get_cost(x):
       params[0].value = x[0]
       copy_parameters(params, "keys", "data")
       apply_parameters(params, "data")

       # Run external solver
       subprocess.run(["python", "run_fedoo_simulation.py"])

       # Compare with DIC fields, reaction forces, etc.
       return compute_cost_from_results()


Gallery Examples
----------------

Two complete examples demonstrate the identification workflow:

- **Hyperelastic identification**: Mooney-Rivlin parameters from Treloar data
  using ``differential_evolution`` and simcoon stress functions
  (see ``examples/analysis/hyperelastic_parameter_identification.py``)

- **Composite identification**: reinforcement properties from effective
  modulus data using Mori-Tanaka and Self-Consistent homogenization
  (see ``examples/heterogeneous/composite_parameter_identification.py``)
