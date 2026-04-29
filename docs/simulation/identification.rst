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
3. **Optimizer**: ``sim.identification()`` wraps
   ``scipy.optimize.differential_evolution``, or call scipy directly
   for more control

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


``identification()`` — Global Optimization
-------------------------------------------

``sim.identification()`` wraps ``scipy.optimize.differential_evolution``
using the bounds from your ``Parameter`` objects. After optimization, the
identified values are written back to each ``Parameter.value``.

.. code-block:: python

   from simcoon.identify import identification
   from simcoon.parameter import Parameter

   params = [
       Parameter(0, bounds=(10000, 200000), key="@Ef",
                 sim_input_files=["Nellipsoids0.dat"]),
       Parameter(1, bounds=(0.01, 0.45), key="@nuf",
                 sim_input_files=["Nellipsoids0.dat"]),
   ]

   result = identification(my_cost_function, params, seed=42, disp=True)
   print(f"E_f = {params[0].value:.0f}, nu_f = {params[1].value:.3f}")

**Arguments:**

- ``cost_fn``: callable ``f(x) -> float`` where ``x`` is a parameter array
- ``parameters``: list of ``Parameter`` (bounds used for search space)
- ``**kwargs``: forwarded to ``differential_evolution`` (``maxiter``,
  ``popsize``, ``tol``, ``seed``, ``polish=True``, ``disp``, etc.)

**Returns:** ``scipy.optimize.OptimizeResult``

For more control (other optimizers, constraints, custom initialization),
call ``scipy.optimize`` directly — the ``Parameter`` objects provide the
bounds and values you need.


``calc_cost()`` — Multi-Level Weighted Cost
--------------------------------------------

``sim.calc_cost()`` computes a weighted cost function with three levels
of weights, designed for multi-test identification:

.. math::

   C = \text{avg}\left(
       W^{\text{test}}_i \cdot
       W^{\text{resp}}_{i,k} \cdot
       W^{\text{pt}}_{i,k,j} \cdot
       (y^{\exp}_{i,k,j} - y^{\num}_{i,k,j})^2
   \right)

Data is organized as a **list of 2-D arrays**, one per test, each of
shape ``(n_points, n_responses)``:

.. code-block:: python

   from simcoon.identify import calc_cost
   import numpy as np

   # Two tensile tests, each with force + displacement columns
   y_exp = [
       np.column_stack([force_exp_1, disp_exp_1]),  # test 1: (N, 2)
       np.column_stack([force_exp_2, disp_exp_2]),  # test 2: (N, 2)
   ]
   y_num = [
       np.column_stack([force_num_1, disp_num_1]),
       np.column_stack([force_num_2, disp_num_2]),
   ]

   # Simple MSE
   cost = calc_cost(y_exp, y_num)

   # NMSE per response (balances force in N vs disp in mm)
   cost = calc_cost(y_exp, y_num, metric='nmse_per_response')

   # Per-test weights (emphasize test 2)
   cost = calc_cost(y_exp, y_num, w_test=np.array([1.0, 3.0]))

Weight levels
^^^^^^^^^^^^^

Three levels of weights are combined multiplicatively:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Level
     - Argument
     - Description
   * - **Test**
     - ``w_test``
     - ``ndarray (n_tests,)`` — weight per experiment/file.
       Use to emphasize certain tests over others.
   * - **Response**
     - ``w_response``
     - ``list of ndarray (n_responses,)`` — weight per response column.
       To balance responses of different magnitudes, use
       ``metric='nmse_per_response'`` instead of manual weights.
   * - **Point**
     - ``w_point``
     - ``list of ndarray (n_points, n_responses)`` — weight per data point.
       Use for heterogeneous confidence or to mask outliers.

Metrics
^^^^^^^

Built-in metrics (numpy only, no extra dependency):

- ``"mse"`` — Mean Squared Error (default)
- ``"nmse"`` — Normalized MSE (divided by variance of all experimental data)
- ``"nmse_per_response"`` — NMSE computed independently per response column,
  then averaged. Each column is divided by its own ``sum(y_exp^2)``,
  balancing responses of different magnitudes (e.g., force in N vs
  displacement in mm). This is the recommended metric for multi-response
  identification.
- ``"rmse"`` — Root Mean Squared Error
- ``"mae"`` — Mean Absolute Error

With ``scikit-learn`` installed (``pip install simcoon[identify]``):

- ``"r2"`` — R-squared score
- ``"mean_squared_error"``, ``"mean_absolute_error"`` — sklearn wrappers
- Any ``sklearn.metrics`` function that accepts ``sample_weight``

If scikit-learn is not installed and an sklearn metric is requested, a
clear error message with install instructions is shown.


Using with External Solvers
----------------------------

The key system works with any simulation tool. For example, with fedoo
(Finite Element Model Updating — FEMU):

.. code-block:: python

   import subprocess
   from simcoon.parameter import Parameter, copy_parameters, apply_parameters
   from simcoon.identify import identification, calc_cost

   params = [
       Parameter(0, bounds=(100e3, 300e3), key="@E",
                 sim_input_files=["material.json"]),
   ]

   def cost(x):
       params[0].value = x[0]
       copy_parameters(params, "keys", "data")
       apply_parameters(params, "data")

       # Run external solver
       subprocess.run(["python", "run_fedoo_simulation.py"])

       # Load results and compare
       y_num = [np.loadtxt("results/reaction_force.txt")]
       y_exp = [np.loadtxt("exp_data/reaction_force.txt")]
       return calc_cost(y_exp, y_num, metric="nmse_per_response")

   result = identification(cost, params, seed=42)


Gallery Examples
----------------

Two complete examples demonstrate the identification workflow:

- **Hyperelastic identification**: Mooney-Rivlin parameters from Treloar data
  using ``differential_evolution`` and simcoon stress functions
  (see ``examples/analysis/hyperelastic_parameter_identification.py``)

- **Composite identification**: reinforcement properties from effective
  modulus data using Mori-Tanaka and Self-Consistent homogenization
  (see ``examples/heterogeneous/composite_parameter_identification.py``)
