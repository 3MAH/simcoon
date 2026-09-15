Use the solver
================================

The Simcoon solver allows you to simulate the mechanical or thermomechanical response of materials under various loading conditions. This page documents the loading path, kept in a ``path.json`` file, and how to drive the solver with it. Since simcoon 2.0 the C++ engine reads no file at all: ``sim.solver.load_path_json`` reads the path in Python into loading objects, and ``sim.solver.solve`` runs them and returns the whole history as numpy arrays. See :doc:`python_solver` to build the same loading objects directly, without any file, and to convert the pre-2.0 ``path.txt`` / ``material.dat`` inputs.

Elastic tensile test
--------------------

Probably the first thing you would like to do with Simcoon is to simulate the mechanical response corresponding to a simple tension test, considering an elastic isotropic material.

We first import *simcoon* (the Python simulation module of simcoon) and *numpy*:

.. code-block:: python

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
              "mode": "linear",
              "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
              "value": [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
              "time": 30.0,
              "ninc": 100,
              "Dn_init": 1.0,
              "Dn_mini": 0.1
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
     - tuple of 3 floats
     - The three Euler angles ``(psi, theta, phi)``, in degrees, giving the material orientation with respect to the reference basis
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

The legacy text loading path
----------------------------

Before 2.0 the loading path was a text file (typically ``path.txt``) in the ``data`` folder. Nothing in simcoon reads it any more: ``scripts/legacy_to_json.py`` turns it into ``path.json`` once. Its structure is kept here for reference:

General structure
^^^^^^^^^^^^^^^^^

.. code-block:: none

    #Initial_temperature
    <T_init>
    #Number_of_blocks
    <nblock>

    #Block
    <block_number>
    #Loading_type
    <type>
    #Control_type(NLGEOM)
    <control_type>
    #Repeat
    <ncycle>
    #Steps
    <nstep>

    <step definitions...>

Block parameters
^^^^^^^^^^^^^^^^

**#Initial_temperature**: The initial temperature of the simulation (in Kelvin).

**#Number_of_blocks**: Total number of loading blocks.

**#Block**: Block number (starting from 1).

**#Loading_type**: Defines the physical problem to solve:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Value
     - Description
   * - 1
     - Mechanical problem
   * - 2
     - Thermomechanical problem (coupled heat equation)

**#Control_type(NLGEOM)**: Defines the kinematic framework and control variables. NLGEOM (non-linear geometry) is activated for Control_type ≥ 2:

.. list-table::
   :header-rows: 1
   :widths: 10 20 70

   * - Value
     - NLGEOM
     - Description
   * - 1
     - No
     - Infinitesimal strains/stress (small deformations)
   * - 2
     - Yes
     - Finite deformation with Lagrangian control (Green-Lagrange strain :math:`\mathbf{E}` / 2nd Piola-Kirchhoff stress :math:`\mathbf{S}`)
   * - 3
     - Yes
     - Finite deformation with logarithmic (true) strain :math:`\boldsymbol{\varepsilon}` / Kirchhoff stress :math:`\boldsymbol{\tau}`
   * - 4
     - Yes
     - Finite deformation with Biot strain :math:`\mathbf{U} - \mathbf{I}` / Biot stress :math:`\mathbf{T}_B = \frac{1}{2}(\mathbf{R}^T\mathbf{P} + \mathbf{P}^T\mathbf{R})`
   * - 5
     - Yes
     - Finite deformation with deformation gradient :math:`\mathbf{F}` control (Eulerian velocity L)
   * - 6
     - Yes
     - Finite deformation with displacement gradient :math:`\nabla\mathbf{u}` control

**#Repeat**: Number of times the block is repeated (for cyclic loading).

**#Steps**: Number of steps within the block.

Step definitions
^^^^^^^^^^^^^^^^

Each step starts with a mode definition:

**#Mode**: Step mode:

- **1**: Linear evolution
- **2**: Sinusoidal evolution
- **3**: Tabular (from a file)

Linear and sinusoidal steps (Mode 1 and 2)
""""""""""""""""""""""""""""""""""""""""""

.. code-block:: none

    #Mode
    1
    #Dn_init 1.
    #Dn_mini 0.1
    #Dn_inc 0.01
    #time
    30.
    #mechanical_state
    E 0.01 
    S 0 S 0
    S 0 S 0 S 0
    #temperature_state
    T 293.5

Parameters:

- **#Dn_init**: Initial size of the first increment (usually 1.0)
- **#Dn_mini**: Minimal size of an increment for convergence issues
- **#Dn_inc**: Increment size as a fraction of the step (0.01 means 100 increments)
- **#time**: Duration of the step :math:`\Delta t`. The time increment is :math:`\delta t = \Delta t \times \delta n`

Mechanical state specification
""""""""""""""""""""""""""""""

For **Control_type = 1** (infinitesimal strains), components are organized in symmetric lower triangular form:

.. code-block:: none

    11
    12 22
    13 23 33

The letter **'S'** indicates stress control, **'E'** indicates strain control:

.. code-block:: none

    E 0.01      # E_11 = 0.01 (strain controlled)
    S 0 S 0     # S_12 = 0, S_22 = 0 (stress controlled)
    S 0 S 0 S 0 # S_13 = 0, S_23 = 0, S_33 = 0 (stress controlled)

For **Control_type = 2, 3, 4** (finite deformation with Lagrangian or logarithmic control), the same symmetric format is used for strain/stress components, with an additional **#spin** block for control types 2, 3, and 4:

.. code-block:: none

    #mechanical_state
    S 3.
    S 0 S 0
    S 0 S 0 S 0
    #spin
    0. 0. 0.
    0. 0. 0.
    0. 0. 0.

The spin tensor :math:`\mathbf{W}` is specified as a full 3×3 matrix.

For **Control_type = 5** (deformation gradient control), the deformation gradient :math:`\mathbf{F}` is specified as a full 3×3 matrix:

.. code-block:: none

    #prescribed_mechanical_state
    5. 0. 0.
    0. 0.4472135955 0.
    0. 0. 0.4472135955

.. note::

   The keywords used as labels (e.g., ``#prescribed_mechanical_state``, ``#prescribed_temperature_state``, ``#mechanical_state``) are placeholders. The solver reads past them and parses the values that follow, so any label can be used.

Temperature state specification
"""""""""""""""""""""""""""""""

For **Loading_type = 1** (mechanical):

.. code-block:: none

    #temperature_state
    T 293.5

The letter **'T'** indicates the temperature at the end of the step.

For **Loading_type = 2** (thermomechanical), additional options are available:

- **T**: Temperature control (imposed temperature)
- **Q**: Heat flux control (imposed heat flux to the RVE)
- **C**: Convection boundary condition

.. code-block:: none

    #prescribed_temperature_state
    Q 0       # Adiabatic conditions (no heat flux)

Tabular steps (Mode 3)
""""""""""""""""""""""

A tabular step follows a table of increments instead of a linear ramp. The
table is the one input that is not JSON: it lives in its own CSV file next to
``path.json`` (``<stem>_tab<k>.csv``, ``k`` numbering the tabular steps of the
path from 1), and the step references it by name:

.. code-block:: json

    {
      "mode": "tabular",
      "tabular": "path_tab1.csv",
      "control": ["strain", "zero", "zero", "zero", "zero", "zero"],
      "Dn_init": 1.0,
      "Dn_mini": 0.01,
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

.. code-block:: none

    #Initial_temperature
    293.15
    #Number_of_blocks
    1

    #Block
    1
    #Loading_type
    1
    #Control_type(NLGEOM)
    1
    #Repeat
    1
    #Steps
    5

    #Mode
    1
    #Dn_init 1.
    #Dn_mini 1.
    #Dn_inc 0.001
    #time
    300
    #prescribed_mechanical_state
    S 1000
    S 0 S 0
    S 0 S 0 S 0
    #prescribed_temperature_state
    T 293.15

    #Mode
    1
    #Dn_init 1.
    #Dn_mini 1.
    #Dn_inc 0.001
    #time
    300
    #prescribed_mechanical_state
    S -1100
    S 0 S 0
    S 0 S 0 S 0
    #prescribed_temperature_state
    T 293.15

    ... (additional steps for cyclic loading)

Hyperelasticity with deformation gradient control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    #Initial_temperature
    293.5
    #Number_of_blocks
    1

    #Block
    1
    #Loading_type
    1
    #Control_type(NLGEOM)
    5
    #Repeat
    1
    #Steps
    1

    #Mode
    1
    #Dn_init 1.
    #Dn_mini 1.
    #Dn_inc 0.1
    #time
    5.
    #prescribed_mechanical_state
    5. 0. 0.
    0. 0.4472135955 0.
    0. 0. 0.4472135955
    #prescribed_temperature_state
    T 290

This applies a uniaxial stretch with :math:`\lambda_1 = 5` and :math:`\lambda_2 = \lambda_3 = 1/\sqrt{5}` (incompressible).

Finite deformation with spin (logarithmic strain)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    #Initial_temperature
    290
    #Number_of_blocks
    1

    #Block
    1
    #Loading_type
    1
    #Control_type(NLGEOM)
    3
    #Repeat
    1
    #Steps
    1

    #Mode
    1
    #Dn_init 1.
    #Dn_mini 1
    #Dn_inc 0.01
    #time
    30.
    #mechanical_state
    S 3.
    S 0 S 0
    S 0 S 0 S 0
    #spin
    0. 0. 0.
    0. 0. 0.
    0. 0. 0.
    #temperature_state
    T 293.5

Thermomechanical loading
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    #Initial_temperature
    290
    #Number_of_blocks
    1

    #Block
    1
    #Loading_type
    2
    #Control_type(NLGEOM)
    1
    #Repeat
    1
    #Steps
    2

    #Mode
    1
    #Dn_init 1.
    #Dn_mini 1.
    #Dn_inc 0.01
    #time
    1
    #prescribed_mechanical_state
    E 0.02
    S 0 S 0
    S 0 S 0 S 0
    #prescribed_temperature_state
    Q 0

    #Mode
    1
    #Dn_init 1.
    #Dn_mini 1
    #Dn_inc 0.01
    #time
    1
    #prescribed_mechanical_state
    E 0.
    S 0 S 0
    S 0 S 0 S 0
    #prescribed_temperature_state
    Q 0

This simulates a strain-controlled loading followed by unloading under adiabatic conditions (Q = 0).