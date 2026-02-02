Solver Architecture Comparison
==============================

Simcoon provides three solver implementations for material point simulations. This document compares their architectures, APIs, and performance characteristics.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
--------

+----------------------+--------------------+--------------------+--------------------+
| Feature              | Legacy C++ Solver  | Python Solver      | C++ Optimized      |
+======================+====================+====================+====================+
| **Location**         | master branch      | feature branch     | feature branch     |
+----------------------+--------------------+--------------------+--------------------+
| **Input Format**     | Text files         | Python objects     | Python objects     |
+----------------------+--------------------+--------------------+--------------------+
| **Output Format**    | Text files         | Python objects     | Python objects     |
+----------------------+--------------------+--------------------+--------------------+
| **UMAT Dispatch**    | Dynamic map        | Via umat_inplace() | Static singleton   |
+----------------------+--------------------+--------------------+--------------------+
| **Buffer Allocation**| Per-increment      | Per-solve          | Pre-allocated      |
+----------------------+--------------------+--------------------+--------------------+
| **Debuggability**    | GDB/LLDB only      | Python debugger    | GDB/LLDB only      |
+----------------------+--------------------+--------------------+--------------------+
| **Extensibility**    | Limited            | Full Python        | Limited            |
+----------------------+--------------------+--------------------+--------------------+


Legacy C++ Solver (Baseline Reference)
--------------------------------------

The legacy solver uses file-based I/O and was the original implementation.

**Input Files Required:**

* ``path.txt`` - Loading path definition
* ``material.dat`` - Material properties
* ``output.dat`` - Output configuration
* ``solver_essentials.inp`` - Solver type and corate
* ``solver_control.inp`` - Newton-Raphson parameters

**API:**

.. code-block:: cpp

   // C++ API
   simcoon::solver(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve,
                   solver_type, corate_type, div_tnew_dt_solver,
                   mul_tnew_dt_solver, miniter_solver, maxiter_solver,
                   inforce_solver, precision_solver, lambda_solver,
                   path_data, path_results, pathfile, outputfile);

**Characteristics:**

* Approximately 1100 lines of C++ code
* Dynamic UMAT name→function mapping (rebuilt per call)
* Newton-Raphson buffers allocated per increment
* File I/O overhead for input/output
* Supports mechanical (type 1) and thermomechanical (type 2) blocks


Python Solver
-------------

The Python solver provides a flexible, object-oriented API with Python control flow.

**API:**

.. code-block:: python

   from simcoon.solver import Solver, Block, StepMeca
   import numpy as np

   # Define material and loading
   props = np.array([210000.0, 0.3, 0.0])  # E, nu, alpha
   step = StepMeca(
       DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
       control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
       Dn_init=100
   )
   block = Block(
       steps=[step],
       umat_name='ELISO',
       props=props,
       nstatev=1
   )

   # Solve
   solver = Solver(blocks=[block], max_iter=10, tol=1e-9)
   history = solver.solve()

   # Access results directly
   final = history[-1]
   print(f"Final stress: {final.sigma}")
   print(f"Final strain: {final.Etot}")

**Characteristics:**

* Approximately 500 lines of Python code
* In-memory data structures (Block, Step, HistoryPoint)
* Full Python debugging support (breakpoints, step-through)
* UMAT calls via ``scc.umat_inplace()`` crossing Python/C++ boundary
* Easy to extend and customize


C++ Optimized Solver
--------------------

The C++ optimized solver provides the same API as the Python solver but with full C++ execution.

**API:**

.. code-block:: python

   from simcoon.solver import Block, StepMeca
   import simcoon._core as scc
   import numpy as np

   # Same setup as Python Solver
   props = np.array([210000.0, 0.3, 0.0])
   step = StepMeca(
       DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
       control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
       Dn_init=100
   )
   block = Block(
       steps=[step],
       umat_name='ELISO',
       props=props,
       nstatev=1
   )

   # Use optimized solver
   history = scc.solver_optimized(blocks=[block], max_iter=10, tol=1e-9)

   # Same HistoryPoint interface
   final = history[-1]
   print(f"Final stress: {final.sigma}")

**Characteristics:**

* Approximately 800 lines of C++ code
* Static UMAT dispatch via singleton (O(1) lookup, built once)
* Pre-allocated Newton-Raphson buffers (reused across increments)
* No Python interpreter overhead during solve loop
* Same output format as Python solver (List[HistoryPoint])


Performance Comparison
----------------------

Benchmark conditions: 100 increments, 10 repetitions, single material point.
All times relative to **Legacy C++ Solver as baseline (1.0x)**.

+----------------------+------------+----------------+------------------+
| UMAT                 | Legacy C++ | Python Solver  | C++ Optimized    |
+======================+============+================+==================+
| **ELISO** (elastic)  | 1.0x       | **7.2x faster**| **46.7x faster** |
|                      | (1.80 ms)  | (0.25 ms)      | (0.04 ms)        |
+----------------------+------------+----------------+------------------+
| **EPICP** (plastic)  | 1.0x       | **2.3x faster**| **4.8x faster**  |
|                      | (1.80 ms)  | (0.80 ms)      | (0.38 ms)        |
+----------------------+------------+----------------+------------------+
| **NEOHC** (finite)   | 1.0x       | **2.9x faster**| **18.6x faster** |
|                      | (1.80 ms)  | (0.62 ms)      | (0.10 ms)        |
+----------------------+------------+----------------+------------------+

**Key Observations:**

* **Legacy C++ solver is dominated by file I/O overhead** (~1.8 ms constant)
* **Both new solvers are faster than legacy** due to in-memory data structures
* **C++ Optimized provides 5-47x speedup** over legacy depending on material
* **Elastoplastic (EPICP) shows smallest gains** because UMAT computation dominates
* **Elastic/hyperelastic show largest gains** because solver overhead dominates


Why Are New Solvers Faster Than Legacy?
---------------------------------------

The legacy C++ solver's ~1.8 ms execution time is dominated by file I/O overhead:

1. **Reading input files** (path.txt, material.dat, output.dat, solver configs)
2. **Writing output files** (results_job_global.txt, results_job_local.txt)
3. **Dynamic UMAT dispatch** (std::map rebuilt on each call)

The new solvers eliminate these bottlenecks:

* **In-memory data structures**: Block/Step/HistoryPoint objects replace file I/O
* **Static UMAT dispatch**: Singleton pattern with O(1) lookup (C++ Optimized)
* **Pre-allocated buffers**: Newton-Raphson arrays reused across increments (C++ Optimized)


Recommendations
---------------

**Use Python Solver when:**

* Developing and debugging new simulations
* Prototyping loading conditions
* Custom post-processing during simulation
* Learning the API
* Need to inspect intermediate states

**Use C++ Optimized Solver when:**

* Running production simulations
* Performance is critical (up to 47x faster than legacy)
* Many simulations in parameter studies
* Integration with optimization loops
* Batch processing of multiple load cases

**Example: Parameter Study**

.. code-block:: python

   import numpy as np
   import simcoon._core as scc
   from simcoon.solver import Block, StepMeca

   # Parameter sweep using optimized solver
   E_values = np.linspace(50000, 200000, 20)
   results = []

   for E in E_values:
       block = Block(
           steps=[StepMeca(DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
                          control=['strain'] + ['stress']*5,
                          Dn_init=100)],
           umat_name='ELISO',
           props=np.array([E, 0.3, 0.0]),
           nstatev=1
       )
       history = scc.solver_optimized(blocks=[block])
       results.append(history[-1].sigma[0])

   print(f"Stress range: {min(results):.2f} - {max(results):.2f}")


Migration from Legacy Solver
----------------------------

To migrate from the legacy file-based solver:

1. **Replace file input with Python objects:**

   .. code-block:: python

      # Legacy (path.txt)
      # #Consigne
      # E 0.02
      # S 0 S 0
      # S 0 S 0 S 0

      # New (Python)
      step = StepMeca(
          DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
          control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress']
      )

2. **Use Python objects instead of file paths:**

   .. code-block:: python

      # Legacy
      scc.solver('ELISO', props, nstatev, 0, 0, 0, 0, 0,
                 'data', 'results', 'path.txt', 'output.txt')

      # New
      history = scc.solver_optimized(blocks=[block])

3. **Access results directly:**

   .. code-block:: python

      # Legacy: Parse output files
      # New: Direct access
      for point in history:
          print(f"Strain: {point.Etot}, Stress: {point.sigma}")
