Use the solver
================================

The Simcoon solver allows you to simulate the mechanical or thermomechanical response of materials under various loading conditions. This page documents the loading path, kept in a ``path.json`` file, and how to drive the solver with it. Since simcoon 2.0 the C++ engine reads no file at all: ``sim.solver.load_path_json`` reads the path in Python into loading objects, and ``sim.solver.solve`` runs them and returns the whole history as numpy arrays. See :doc:`python_solver` to build the same loading objects directly, without any file, and to convert the pre-2.0 ``path.txt`` / ``material.dat`` inputs.

Elastic tensile test
--------------------

Probably the first thing you would like to do with Simcoon is to simulate the mechanical response corresponding to a simple tension test, considering an elastic isotropic material.

We first import *simcoon* (the Python simulation module of simcoon) and *numpy*:

.. code-block:: python

    import os

    import numpy as np
    import simcoon as sim

Next we shall define the material constitutive law to be utilized and the associated material properties. We will pass them as a numpy array:

.. code-block:: python

    umat_name = 'ELISO'  # This is the 5 character code for the elastic-isotropic subroutine
    nstatev = 1  # The number of scalar state variables required

    E = 700000.  # The Young modulus
    nu = 0.2  # The Poisson coefficient
    alpha = 1.E-5  # The coefficient of thermal expansion

    # Three Euler angles to represent the material orientation with respect to the reference basis
    psi_rve = 0.
    theta_rve = 0.
    phi_rve = 0.

    # Solver parameters
    solver_type = 0  # Solver strategy (0 for Newton-Raphson)
    corate_type = 3  # Corotational spin rate type (0: Jaumann, 1: Green-Naghdi, 2: logarithmic/XBM, 3: logarithmic_R — recommended default)

    props = np.array([E, nu, alpha])

We shall then define the location of the loading path file:

.. code-block:: python

    path_data = 'data'
    pathfile = 'path.json'

The last part is to define the loading path. Create a folder ``data`` and a file named ``path.json`` with the following content:

.. code-block:: json

    {
      "initial_temperature": 293.5,
      "corate": "logarithmic_R",
      "blocks": [
        {
          "control_type": "small_strain",
          "ncycle": 1,
          "steps": [
            {
              "thermomechanical": false,
              "mode": "linear",
              "time": 30.0,
              "ninc": 100,
              "Dn_init": 1.0,
              "Dn_mini": 0.1,
              "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
              "value": [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
              "T_final": 293.5
            }
          ]
        }
      ]
    }

This corresponds to a pure strain-controlled tension test in direction 1 up to 1% strain, at 293.5K, in 100 increments. The file is what ``sim.solver.save_path_json`` writes for a :class:`~simcoon.solver.Block` of one :class:`~simcoon.solver.StepMeca`; the pre-2.0 text format is converted to it once with ``scripts/legacy_to_json.py`` (see :doc:`python_solver`).

Finally, read the path file and run the solver:

.. code-block:: python

    blocks, T_init, _ = sim.solver.load_path_json(os.path.join(path_data, pathfile))

    res = sim.solver.solve(
        blocks,
        umat_name,
        props,
        nstatev,
        T_init=T_init,
        solver_type=solver_type,
        corate=corate_type,
        orientation=(psi_rve, theta_rve, phi_rve),
    )

No result file is written: ``res`` holds the whole history in memory, as numpy
arrays of shape ``(6, N)`` for the tensor quantities:

.. code-block:: python

    e11, e22, e33, e12, e13, e23 = res["Strain"]   # total strain
    s11, s22, s33, s12, s13, s23 = res["Stress"]   # Cauchy stress
    time, T = res["Time"], res["Temp"]
    Wm, Wm_r, Wm_ir, Wm_d = res["Wm"]

See :doc:`python_solver` for the full list of available keys.

Solver parameters
-----------------

``sim.solver.solve`` takes the following parameters:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Type
     - Description
   * - blocks
     - Block, StepMeca, or a sequence of them
     - The loading path, as returned by ``sim.solver.load_path_json`` or built directly (see :doc:`python_solver`)
   * - umat_name
     - string or callable
     - 5-character code identifying the constitutive law (e.g., 'ELISO', 'EPICP', 'EPKCP'), or a constitutive law written in Python
   * - props
     - numpy array
     - Material properties array
   * - nstatev
     - int
     - Number of internal state variables
   * - T_init
     - float
     - Initial temperature in Kelvin (default 293.15); ``load_path_json`` returns the value read from the path file
   * - corate
     - int or string
     - Corotational spin rate type (see below)
   * - tangent_mode
     - int
     - Tangent-operator mode (default 2): see below
   * - solver_type
     - int
     - Solver strategy (0: Newton-Raphson)
   * - orientation
     - simcoon.Rotation or 3 floats
     - Orientation of the material frame: a ``Rotation`` or its Euler angles ``(psi, theta, phi)`` in degrees (see :doc:`python_solver`, mean-field composites, for the convention)
   * - phases
     - sequence, optional
     - Sub-phases of a mean-field model (MIHEN, MIMTN, MISCN, MIPLN), passed in memory: :class:`~simcoon.solver.micromechanics.Ellipsoid` / :class:`~simcoon.solver.micromechanics.Layer` objects (numbered by position) or the dicts :func:`~simcoon.solver.micromechanics.to_phase_dicts` makes of them. Their geometry must be the one the model builds, the concentrations must sum to 1, and a sub-phase that is itself a mean-field model carries its own sub-phases in its ``phases`` attribute
   * - record_tangent
     - bool
     - Whether to record the tangent operator in the results (default True)

Tangent-operator modes
^^^^^^^^^^^^^^^^^^^^^^

The ``tangent_mode`` keyword selects the tangent operator returned by the
constitutive models (also exposed as named constants:
``sim.tangent_none``, ``sim.tangent_continuum``, ``sim.tangent_algorithmic``,
``sim.tangent_closest_point``):

.. list-table::
   :header-rows: 1
   :widths: 10 25 65

   * - Value
     - Name
     - Description
   * - 0
     - none
     - No tangent assembly — ``Lt`` stays the elastic operator (explicit
       integration schemes)
   * - 1
     - continuum
     - Continuum elasto-(visco)plastic operator (this was mode 0 before the
       2.0 renumbering)
   * - 2
     - algorithmic
     - Simo–Hughes consistent (algorithmic) operator — **default**; exact
       Jacobian of the discrete return map for J2-type flows, Q-quadratic
       global convergence (this was mode 1 before the 2.0 renumbering)
   * - 3
     - closest-point
     - Reserved for the closest-point-projection exact operator (future
       release); currently raises an error

.. note::
   **2.0 renumbering.** Pre-2.0, ``tangent_mode 0`` meant *continuum* and
   ``1`` meant *algorithmic*; there was no "none" mode. The converged
   response is identical in every mode (the tangent steers the global
   Newton iteration, not the residual) — only iteration counts and run
   time change. Scripts passing explicit values should shift them by +1;
   scripts relying on the default silently upgrade from continuum to the
   (faster) algorithmic operator.

Corotational spin rate types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``corate`` parameter controls the corotational formulation used in finite deformation problems:

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
     - Uses the logarithmic (Xi--Meyers--Brühns) spin rate
   * - 3
     - Logarithmic_R (log_R)
     - Logarithmic strain transported by the exact polar rotation increment
       :math:`\Delta\mathbf{R} = \mathbf{R}_1\mathbf{R}_0^T` (recommended
       default: exact tangent transport, including with rotated
       internal-variable history)
   * - 4
     - Truesdell
     - Convected (Truesdell) rate, frame increment :math:`\Delta\mathbf{F}`
   * - 5
     - Logarithmic_F (log_F)
     - Convected logarithmic rate (pure :math:`\mathbf{F}` transport)

The loading path file
---------------------

``path.json`` is what ``sim.solver.save_path_json`` writes for a list of
:class:`~simcoon.solver.Block` objects, and what ``sim.solver.load_path_json``
reads back. Every name below is a keyword of :class:`~simcoon.solver.Block`,
:class:`~simcoon.solver.StepMeca` or :class:`~simcoon.solver.StepThermomeca`, so
a path built in Python and a path read from the file are the same objects (see
:doc:`python_solver`). Enumerated entries (``control_type``, ``mode``,
``control``, ``thermal_control``, ``corate``) accept the names given here or the
integer codes of the C++ solver; ``save_path_json`` always writes the names.

General structure
^^^^^^^^^^^^^^^^^

.. code-block:: json

    {
      "initial_temperature": 293.15,
      "corate": "logarithmic_R",
      "blocks": [
        {
          "control_type": "small_strain",
          "ncycle": 1,
          "steps": [
            { "...": "step definitions" }
          ]
        }
      ]
    }

Path and block parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

**initial_temperature**: initial temperature of the simulation (Kelvin); ``load_path_json`` returns it as ``T_init`` for ``solve``.

**corate**: objective rate of the finite-strain control types (``"jaumann"``, ``"green_naghdi"``, ``"logarithmic"``, ``"logarithmic_R"``, ``"truesdell"``, ``"logarithmic_F"``, see the table above); returned by ``load_path_json`` for ``solve(corate=...)``.

**blocks**: the loading blocks, run in order. A block is mechanical or thermomechanical (coupled heat equation) according to its steps: all the steps of a block are :class:`~simcoon.solver.StepMeca` (``"thermomechanical": false``) or all are :class:`~simcoon.solver.StepThermomeca` (``"thermomechanical": true``).

**control_type**: kinematic framework and control variables of the block. NLGEOM (non-linear geometry) is activated from ``"green_lagrange"`` on; thermomechanical blocks only support ``"small_strain"``:

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Name (code)
     - NLGEOM
     - Description
   * - ``"small_strain"`` (1)
     - No
     - Infinitesimal strains/stress (small deformations)
   * - ``"green_lagrange"`` (2)
     - Yes
     - Finite deformation with Lagrangian control (Green-Lagrange strain :math:`\mathbf{E}` / 2nd Piola-Kirchhoff stress :math:`\mathbf{S}`)
   * - ``"logarithmic"`` (3)
     - Yes
     - Finite deformation with logarithmic (true) strain :math:`\boldsymbol{\varepsilon}` / Kirchhoff stress :math:`\boldsymbol{\tau}`
   * - ``"biot"`` (4)
     - Yes
     - Finite deformation with the right stretch :math:`\mathbf{U}` / Biot stress :math:`\mathbf{T}_B = \frac{1}{2}(\mathbf{R}^T\mathbf{P} + \mathbf{P}^T\mathbf{R})`. The kinematic targets are the components of :math:`\mathbf{U}` itself, not of the Biot strain :math:`\mathbf{U} - \mathbf{I}`: the undeformed state is 1 on the diagonal, and a 8 % stretch is ``1.08``
   * - ``"F"`` (5)
     - Yes
     - Finite deformation with deformation gradient :math:`\mathbf{F}` control (Eulerian velocity L)
   * - ``"gradU"`` (6)
     - Yes
     - Finite deformation with displacement gradient :math:`\nabla\mathbf{u}` control

**ncycle**: number of times the step sequence of the block is repeated (cyclic loading). A block holding a tabular step cannot be cycled.

**steps**: the steps of the block, in order.

Step definitions
^^^^^^^^^^^^^^^^

**mode**: evolution of the prescribed components over the step:

- ``"linear"`` (1): linear ramp to the target values
- ``"sinusoidal"`` (2): sinusoidal evolution to the target values
- ``"tabular"`` (3): table of increments read from a CSV file

Linear and sinusoidal steps
"""""""""""""""""""""""""""

.. code-block:: json

    {
      "thermomechanical": false,
      "mode": "linear",
      "time": 30.0,
      "ninc": 100,
      "Dn_init": 1.0,
      "Dn_mini": 0.1,
      "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
      "value": [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
      "T_final": 293.5
    }

Parameters:

- **time**: duration of the step :math:`\Delta t`
- **ninc**: number of increments; the time increment is :math:`\delta t = \Delta t / n_{inc}`
- **Dn_init**: initial size of the first sub-increment, as a fraction of one increment (usually 1.0)
- **Dn_mini**: minimal sub-increment fraction the adaptive stepping may cut down to before giving up
- **control**, **value**: the prescribed mechanical state (below)
- **T_final**: temperature at the end of the step, ``null`` to hold the current one (mechanical steps: imposed temperature, thermal expansion without the heat equation)

Mechanical state specification
""""""""""""""""""""""""""""""

``control`` names the prescribed quantity of each component and ``value`` its target at the end of the step, in Voigt order:

.. code-block:: none

    11  22  33  12  13  23

``"strain"`` (or ``"E"``) prescribes the kinematic quantity of the control type (strain component), ``"stress"`` (or ``"S"``) the static one (stress component); a single string applies to all six components. For example ``"control": ["strain", "stress", "stress", "stress", "stress", "stress"]`` with ``"value": [0.01, 0, 0, 0, 0, 0]`` is a uniaxial tension test to 1 % strain in direction 1, the other components being stress-free. Targets are absolute values, reached at the end of the step whatever the state at its start.

For ``"F"`` and ``"gradU"`` control types the state is the full tensor, 9 components row-major:

.. code-block:: none

    11  12  13  21  22  23  31  32  33

``control`` is then ``"strain"`` for all of them and ``value`` holds the target :math:`\mathbf{F}` (or :math:`\nabla\mathbf{u}`).

For the mixed finite-strain control types (``"green_lagrange"``, ``"logarithmic"``, ``"biot"``) a rotation rate can be superimposed with **BC_w**, a 3x3 spin matrix applied during the step:

.. code-block:: json

    "BC_w": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

Thermal state specification
"""""""""""""""""""""""""""

Mechanical steps only carry **T_final** (above). Thermomechanical steps (``"thermomechanical": true``) solve the heat equation and take **thermal_control**:

- ``"temperature"`` (T): imposed temperature, ramped to **T_final**
- ``"heat_flux"`` (Q): imposed heat flux **Q** on the RVE (``0.0`` for adiabatic conditions)
- ``"convection"`` (C): 0D convection :math:`Q = -q_{conv}\,(T - T_{init})` with coefficient **q_conv**

.. code-block:: json

    {
      "thermomechanical": true,
      "mode": "linear",
      "time": 1.0,
      "ninc": 100,
      "Dn_init": 1.0,
      "Dn_mini": 1.0,
      "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
      "value": [0.02, 0.0, 0.0, 0.0, 0.0, 0.0],
      "thermal_control": "heat_flux",
      "Q": 0.0
    }

Tabular steps
"""""""""""""

A tabular step follows a table of increments instead of a linear ramp. The
table is the one input that is not JSON: it lives in its own CSV file next to
``path.json`` (``<stem>_tab<k>.csv``, ``k`` numbering the tabular steps of the
path from 1), and the step references it by name:

.. code-block:: json

    {
      "thermomechanical": false,
      "mode": "tabular",
      "tabular": "path_tab1.csv",
      "Dn_init": 1.0,
      "Dn_mini": 0.01,
      "control": ["strain", "zero", "zero", "zero", "zero", "zero"],
      "tabular_T": false
    }

The ``control`` list says which components the table drives:

- ``"strain"``: strain-controlled component (a column of the table)
- ``"stress"``: stress-controlled component (a column of the table)
- ``"zero"``: component held at zero (no column)

``tabular_T``: ``false`` if the temperature is constant, ``true`` if the table
carries a temperature column.

The table, one row per increment, comma- or whitespace-separated, ``#`` lines
ignored:

.. code-block:: none

    # time, E11
    0.01, 0.0005
    0.02, 0.0010
    0.03, 0.0015
    ...

Columns: **time** (absolute simulation time, continuing from the previous step),
**T** if ``tabular_T`` is set (**Q** for a heat-flux thermomechanical step),
then the controlled components in Voigt order 11, 22, 33, 12, 13, 23.
``sim.solver.save_path_json`` writes this file (with the header) from the
``tabular`` array of a :class:`~simcoon.solver.StepMeca`; the pre-2.0 ``#File``
increment tables (leading increment number, time restarting at 0) are converted
to it by ``scripts/legacy_to_json.py``. A tabular step cannot be cycled.

Examples
--------

Cyclic loading (plasticity)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Stress-controlled tension/compression cycle, 1000 increments per step:

.. code-block:: json

    {
      "initial_temperature": 293.15,
      "corate": "logarithmic_R",
      "blocks": [
        {
          "control_type": "small_strain",
          "ncycle": 1,
          "steps": [
            {
              "thermomechanical": false,
              "mode": "linear",
              "time": 300.0,
              "ninc": 1000,
              "Dn_init": 1.0,
              "Dn_mini": 1.0,
              "control": ["stress", "stress", "stress", "stress", "stress", "stress"],
              "value": [1000.0, 0.0, 0.0, 0.0, 0.0, 0.0],
              "T_final": 293.15
            },
            {
              "thermomechanical": false,
              "mode": "linear",
              "time": 300.0,
              "ninc": 1000,
              "Dn_init": 1.0,
              "Dn_mini": 1.0,
              "control": ["stress", "stress", "stress", "stress", "stress", "stress"],
              "value": [-1100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
              "T_final": 293.15
            }
          ]
        }
      ]
    }

(additional steps, or ``"ncycle": 10`` on the block, for cyclic loading; the
files are written exactly like this by ``sim.solver.save_path_json``)

Hyperelasticity with deformation gradient control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

    {
      "initial_temperature": 293.5,
      "corate": "logarithmic_R",
      "blocks": [
        {
          "control_type": "F",
          "ncycle": 1,
          "steps": [
            {
              "thermomechanical": false,
              "mode": "linear",
              "time": 5.0,
              "ninc": 10,
              "Dn_init": 1.0,
              "Dn_mini": 1.0,
              "control": ["strain", "strain", "strain",
                          "strain", "strain", "strain",
                          "strain", "strain", "strain"],
              "value": [5.0, 0.0, 0.0,
                        0.0, 0.4472135955, 0.0,
                        0.0, 0.0, 0.4472135955],
              "T_final": 290.0
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
      "initial_temperature": 290.0,
      "corate": "logarithmic_R",
      "blocks": [
        {
          "control_type": "logarithmic",
          "ncycle": 1,
          "steps": [
            {
              "thermomechanical": false,
              "mode": "linear",
              "time": 30.0,
              "ninc": 100,
              "Dn_init": 1.0,
              "Dn_mini": 1.0,
              "control": ["stress", "stress", "stress", "stress", "stress", "stress"],
              "value": [3.0, 0.0, 0.0, 0.0, 0.0, 0.0],
              "BC_w": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
              "T_final": 293.5
            }
          ]
        }
      ]
    }

Thermomechanical loading
^^^^^^^^^^^^^^^^^^^^^^^^

Adiabatic strain-controlled loading/unloading (heat equation on, no heat flux):

.. code-block:: json

    {
      "initial_temperature": 290.0,
      "corate": "logarithmic_R",
      "blocks": [
        {
          "control_type": "small_strain",
          "ncycle": 1,
          "steps": [
            {
              "thermomechanical": true,
              "mode": "linear",
              "time": 1.0,
              "ninc": 100,
              "Dn_init": 1.0,
              "Dn_mini": 1.0,
              "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
              "value": [0.02, 0.0, 0.0, 0.0, 0.0, 0.0],
              "thermal_control": "heat_flux",
              "Q": 0.0
            },
            {
              "thermomechanical": true,
              "mode": "linear",
              "time": 1.0,
              "ninc": 100,
              "Dn_init": 1.0,
              "Dn_mini": 1.0,
              "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
              "value": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
              "thermal_control": "heat_flux",
              "Q": 0.0
            }
          ]
        }
      ]
    }

Pre-2.0 text inputs
^^^^^^^^^^^^^^^^^^^

Before 2.0 the loading path was a text file (``path.txt``, with ``#Mode`` /
``#prescribed_mechanical_state`` blocks and ``#File`` increment tables). Nothing
in simcoon reads that format any more: ``scripts/legacy_to_json.py`` converts a
``data`` directory once into the files described here (see :doc:`python_solver`).
