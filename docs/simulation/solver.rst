Use the solver
================================

The Simcoon solver allows you to simulate the mechanical or thermomechanical response of materials under various loading conditions. This documentation covers the Python interface, solver parameters, and JSON-based configuration.

Elastic tensile test
--------------------

Probably the first thing you would like to do with Simcoon is to simulate the mechanical response corresponding to a simple tension test, considering an elastic isotropic material.

We first import *simcoon* (the Python simulation module of simcoon) and *numpy*:

.. code-block:: python

    import numpy as np
    import simcoon as sim
    from simcoon.solver.io import load_material_json, load_path_json

Material Configuration (JSON)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a file ``data/material.json``:

.. code-block:: json

   {
     "name": "ELISO",
     "props": {"E": 700000, "nu": 0.2, "alpha": 1e-5},
     "nstatev": 1,
     "orientation": {"psi": 0, "theta": 0, "phi": 0}
   }

Path Configuration (JSON)
^^^^^^^^^^^^^^^^^^^^^^^^^

Create a file ``data/path.json`` for a simple tension test up to 1% strain:

.. code-block:: json

   {
     "initial_temperature": 293.5,
     "blocks": [
       {
         "type": "mechanical",
         "control_type": "small_strain",
         "corate_type": "logarithmic",
         "ncycle": 1,
         "steps": [
           {
             "time": 30.0,
             "Dn_init": 1.0,
             "Dn_mini": 0.1,
             "Dn_inc": 0.01,
             "DEtot": [0.01, 0, 0, 0, 0, 0],
             "Dsigma": [0, 0, 0, 0, 0, 0],
             "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
             "DT": 0
           }
         ]
       }
     ]
   }

Running the Simulation
^^^^^^^^^^^^^^^^^^^^^^

**Option 1: Using Python Solver API (Recommended)**

The Python Solver API provides full control over the simulation:

.. code-block:: python

    from simcoon.solver import Solver, Block, StepMeca

    # Material properties: E, nu, alpha
    props = np.array([700000, 0.2, 1e-5])

    # Define loading step
    step = StepMeca(
        DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=1,
        Dn_inc=100,
        time=30.0
    )

    # Create block
    block = Block(
        steps=[step],
        umat_name='ELISO',
        props=props,
        nstatev=1
    )

    # Run solver
    solver = Solver(blocks=[block])
    history = solver.solve()

    # Access results
    for state in history:
        print(f"Strain: {state.Etot[0]:.4f}, Stress: {state.sigma[0]:.2f}")

**Option 2: Using JSON Configuration**

Load configurations from JSON files and convert to Python objects:

.. code-block:: python

    from simcoon.solver.io import load_material_json, load_path_json

    # Load from JSON
    material = load_material_json('data/material.json')
    path = load_path_json('data/path.json')

    # Convert to Block/Step objects and run
    # (See migration guide for full conversion examples)

Python Solver API
-----------------

The ``Solver`` class is the main entry point for running simulations:

.. code-block:: python

    from simcoon.solver import Solver, Block, StepMeca

    solver = Solver(blocks=[block], T_init=293.15)
    history = solver.solve()

**Solver Parameters:**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Type
     - Description
   * - blocks
     - List[Block]
     - List of loading blocks to execute
   * - T_init
     - float
     - Initial temperature (default: 293.15 K)

**Block Parameters:**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Type
     - Description
   * - steps
     - List[StepMeca]
     - List of loading steps in the block
   * - umat_name
     - str
     - 5-character UMAT code (e.g., 'ELISO', 'EPICP', 'EPKCP')
   * - props
     - np.ndarray
     - Material properties array
   * - nstatev
     - int
     - Number of internal state variables
   * - psi
     - float
     - First Euler angle (degrees) for material orientation (default: 0)
   * - theta
     - float
     - Second Euler angle (degrees) for material orientation (default: 0)
   * - phi
     - float
     - Third Euler angle (degrees) for material orientation (default: 0)
   * - ncycles
     - int
     - Number of times to repeat this block (default: 1)
   * - solver_type
     - int
     - Solver strategy: 0=Newton-Raphson (default: 0)
   * - control_type
     - int
     - Kinematic framework (see Control Types below)
   * - corate_type
     - int
     - Corotational spin rate type (see Corotational Types below)

**Control Types (control_type):**

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Value
     - Name
     - Description
   * - 1
     - small_strain
     - Infinitesimal strain (small deformations)
   * - 2
     - logarithmic
     - Logarithmic (true) strain / Kirchhoff stress
   * - 3
     - deformation_gradient
     - Deformation gradient :math:`\mathbf{F}` control

**StepMeca Parameters:**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Type
     - Description
   * - DEtot_end
     - np.ndarray
     - Target strain increment [E11, E12, E22, E13, E23, E33]
   * - Dsigma_end
     - np.ndarray
     - Target stress increment (optional)
   * - control
     - List[str]
     - Control type per component: 'strain' or 'stress'
   * - time
     - float
     - Step duration (default: 1.0)
   * - Dn_init
     - int
     - Initial increment count (default: 1)
   * - Dn_mini
     - int
     - Minimum increment for convergence (default: 1)
   * - Dn_inc
     - int
     - Total number of increments (default: 100)

Corotational spin rate types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``corate_type`` parameter controls the corotational formulation used in finite deformation problems:

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Value
     - Type
     - Description
   * - 0
     - Jaumann
     - Uses the spin tensor :math:`\mathbf{W}` from the velocity gradient
   * - 1
     - Green-Naghdi
     - Uses the spin from the polar decomposition :math:`\dot{\mathbf{R}}\mathbf{R}^T`
   * - 2
     - Logarithmic
     - Uses the logarithmic spin rate (recommended for large deformations)

JSON Path Configuration
-----------------------

The loading path is defined in a JSON file (e.g., ``path.json``) located in the ``data`` folder.

General Structure
^^^^^^^^^^^^^^^^^

.. code-block:: json

   {
     "initial_temperature": 293.15,
     "blocks": [
       {
         "type": "mechanical",
         "control_type": "small_strain",
         "corate_type": "logarithmic",
         "ncycle": 1,
         "steps": [...]
       }
     ]
   }

Block Parameters
^^^^^^^^^^^^^^^^

**type**: Defines the physical problem to solve:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Value
     - Description
   * - ``"mechanical"``
     - Mechanical problem
   * - ``"thermomechanical"``
     - Thermomechanical problem (coupled heat equation)

**control_type**: Defines the kinematic framework and control variables:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Value
     - Description
   * - ``"small_strain"``
     - Infinitesimal strains/stress (small deformations)
   * - ``"lagrangian"``
     - Finite deformation with Lagrangian control (Green-Lagrange strain :math:`\mathbf{E}` / 2nd Piola-Kirchhoff stress :math:`\mathbf{S}`)
   * - ``"logarithmic"``
     - Finite deformation with logarithmic (true) strain :math:`\boldsymbol{\varepsilon}` / Kirchhoff stress :math:`\boldsymbol{\tau}`
   * - ``"biot"``
     - Finite deformation with Biot strain :math:`\mathbf{U} - \mathbf{I}` / Biot stress
   * - ``"deformation_gradient"``
     - Finite deformation with deformation gradient :math:`\mathbf{F}` control
   * - ``"displacement_gradient"``
     - Finite deformation with displacement gradient :math:`\nabla\mathbf{u}` control

**corate_type**: Corotational spin rate type (``"jaumann"``, ``"green_naghdi"``, or ``"logarithmic"``).

**ncycle**: Number of times the block is repeated (for cyclic loading).

Step Definitions
^^^^^^^^^^^^^^^^

Each step in the ``"steps"`` array contains:

.. code-block:: json

   {
     "time": 30.0,
     "Dn_init": 1.0,
     "Dn_mini": 0.1,
     "Dn_inc": 0.01,
     "DEtot": [0.01, 0, 0, 0, 0, 0],
     "Dsigma": [0, 0, 0, 0, 0, 0],
     "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
     "DT": 0
   }

Step Parameters:

- **time**: Duration of the step :math:`\Delta t`
- **Dn_init**: Initial size of the first increment (usually 1.0)
- **Dn_mini**: Minimal size of an increment for convergence issues
- **Dn_inc**: Increment size as a fraction of the step (0.01 means 100 increments)
- **DEtot**: Strain increments [E11, E12, E22, E13, E23, E33]
- **Dsigma**: Stress increments [S11, S12, S22, S13, S23, S33]
- **control**: Control type for each component (``"strain"`` or ``"stress"``)
- **DT**: Temperature increment

Finite Deformation with Spin
""""""""""""""""""""""""""""

For finite deformation control types that require spin specification:

.. code-block:: json

   {
     "time": 30.0,
     "Dn_init": 1.0,
     "Dn_mini": 1.0,
     "Dn_inc": 0.01,
     "DEtot": [0, 0, 0, 0, 0, 0],
     "Dsigma": [3.0, 0, 0, 0, 0, 0],
     "control": ["stress", "stress", "stress", "stress", "stress", "stress"],
     "spin": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
     "DT": 0
   }

Deformation Gradient Control
""""""""""""""""""""""""""""

For deformation gradient control:

.. code-block:: json

   {
     "time": 5.0,
     "Dn_init": 1.0,
     "Dn_mini": 1.0,
     "Dn_inc": 0.1,
     "F": [[5.0, 0, 0], [0, 0.4472135955, 0], [0, 0, 0.4472135955]],
     "DT": 0
   }

Temperature Control
"""""""""""""""""""

For thermomechanical loading:

- **DT**: Temperature increment for temperature control
- **DQ**: Heat flux for adiabatic conditions (set to 0 for adiabatic)

.. code-block:: json

   {
     "time": 1.0,
     "Dn_init": 1.0,
     "Dn_mini": 1.0,
     "Dn_inc": 0.01,
     "DEtot": [0.02, 0, 0, 0, 0, 0],
     "Dsigma": [0, 0, 0, 0, 0, 0],
     "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
     "DQ": 0
   }

Examples
--------

Cyclic loading (plasticity)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

   {
     "initial_temperature": 293.15,
     "blocks": [
       {
         "type": "mechanical",
         "control_type": "small_strain",
         "corate_type": "jaumann",
         "ncycle": 1,
         "steps": [
           {
             "time": 300,
             "Dn_init": 1.0,
             "Dn_mini": 1.0,
             "Dn_inc": 0.001,
             "DEtot": [0, 0, 0, 0, 0, 0],
             "Dsigma": [1000, 0, 0, 0, 0, 0],
             "control": ["stress", "stress", "stress", "stress", "stress", "stress"],
             "DT": 0
           },
           {
             "time": 300,
             "Dn_init": 1.0,
             "Dn_mini": 1.0,
             "Dn_inc": 0.001,
             "DEtot": [0, 0, 0, 0, 0, 0],
             "Dsigma": [-1100, 0, 0, 0, 0, 0],
             "control": ["stress", "stress", "stress", "stress", "stress", "stress"],
             "DT": 0
           }
         ]
       }
     ]
   }

Hyperelasticity with deformation gradient control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

   {
     "initial_temperature": 293.5,
     "blocks": [
       {
         "type": "mechanical",
         "control_type": "deformation_gradient",
         "ncycle": 1,
         "steps": [
           {
             "time": 5.0,
             "Dn_init": 1.0,
             "Dn_mini": 1.0,
             "Dn_inc": 0.1,
             "F": [[5.0, 0, 0], [0, 0.4472135955, 0], [0, 0, 0.4472135955]],
             "DT": 0
           }
         ]
       }
     ]
   }

This applies a uniaxial stretch with :math:`\lambda_1 = 5` and :math:`\lambda_2 = \lambda_3 = 1/\sqrt{5}` (incompressible).

Finite deformation with spin (logarithmic strain)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

   {
     "initial_temperature": 290,
     "blocks": [
       {
         "type": "mechanical",
         "control_type": "logarithmic",
         "corate_type": "logarithmic",
         "ncycle": 1,
         "steps": [
           {
             "time": 30.0,
             "Dn_init": 1.0,
             "Dn_mini": 1.0,
             "Dn_inc": 0.01,
             "DEtot": [0, 0, 0, 0, 0, 0],
             "Dsigma": [3.0, 0, 0, 0, 0, 0],
             "control": ["stress", "stress", "stress", "stress", "stress", "stress"],
             "spin": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             "DT": 0
           }
         ]
       }
     ]
   }

Thermomechanical loading
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

   {
     "initial_temperature": 290,
     "blocks": [
       {
         "type": "thermomechanical",
         "control_type": "small_strain",
         "ncycle": 1,
         "steps": [
           {
             "time": 1.0,
             "Dn_init": 1.0,
             "Dn_mini": 1.0,
             "Dn_inc": 0.01,
             "DEtot": [0.02, 0, 0, 0, 0, 0],
             "Dsigma": [0, 0, 0, 0, 0, 0],
             "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
             "DQ": 0
           },
           {
             "time": 1.0,
             "Dn_init": 1.0,
             "Dn_mini": 1.0,
             "Dn_inc": 0.01,
             "DEtot": [0, 0, 0, 0, 0, 0],
             "Dsigma": [0, 0, 0, 0, 0, 0],
             "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
             "DQ": 0
           }
         ]
       }
     ]
   }

This simulates a strain-controlled loading followed by unloading under adiabatic conditions (DQ = 0).

Material JSON Format
--------------------

The ``simcoon.solver.io`` module provides JSON I/O for material configurations:

.. code-block:: json

   {
     "name": "ELISO",
     "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5},
     "nstatev": 1,
     "orientation": {"psi": 0, "theta": 0, "phi": 0}
   }

Usage
^^^^^

.. code-block:: python

   from simcoon.solver.io import (
       load_material_json, save_material_json,
       load_path_json, save_path_json,
   )

   # Load from JSON
   material = load_material_json('material.json')
   path = load_path_json('path.json')

   # Save configurations
   save_material_json('material.json', material)
   save_path_json('path.json', path)

See Also
--------

- :doc:`micromechanics` - Python micromechanics I/O (phases, ellipsoids, layers)
- :doc:`output` - Output file configuration