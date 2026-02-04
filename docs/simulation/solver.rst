Use the solver
================================

The Simcoon solver allows you to simulate the mechanical or thermomechanical response of materials under various loading conditions. This documentation covers the Python interface, solver parameters, and the structure of input files.

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
    corate_type = 2  # Corotational spin rate type (0: Jaumann, 1: Green-Naghdi, 2: logarithmic)

    props = np.array([E, nu, alpha])

We shall then define the location of the data input files and the results output file:

.. code-block:: python

    path_data = 'data'
    path_results = 'results'
    pathfile = 'path.txt'
    outputfile = 'results_ELISO.txt'

The last part is to define the loading path. Create a folder ``data`` and a text file named ``path.txt`` with the following content:

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
    1    
    #Repeat
    1
    #Steps
    1

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

This corresponds to a pure strain-controlled tension test in direction 1 up to 1% strain, at 293.5K.

Finally, call the solver function:

.. code-block:: python

    sim.solver(
        umat_name,
        props,
        nstatev,
        psi_rve,
        theta_rve,
        phi_rve,
        solver_type,
        corate_type,
        path_data,
        path_results,
        pathfile,
        outputfile,
    )

The result file ``results_ELISO.txt`` will be created in the ``results`` folder.

Solver parameters
-----------------

The solver function takes the following parameters:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Type
     - Description
   * - umat_name
     - string
     - 5-character code identifying the constitutive law (e.g., 'ELISO', 'EPICP', 'EPKCP')
   * - props
     - numpy array
     - Material properties array
   * - nstatev
     - int
     - Number of internal state variables
   * - psi_rve
     - float
     - First Euler angle (in degrees) for material orientation
   * - theta_rve
     - float
     - Second Euler angle (in degrees) for material orientation
   * - phi_rve
     - float
     - Third Euler angle (in degrees) for material orientation
   * - solver_type
     - int
     - Solver strategy (0: Newton-Raphson)
   * - corate_type
     - int
     - Corotational spin rate type (see below)
   * - path_data
     - string
     - Path to the folder containing input files
   * - path_results
     - string
     - Path to the folder for output files
   * - pathfile
     - string
     - Name of the loading path file (default: 'path.txt')
   * - outputfile
     - string
     - Name of the output result file

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

Define the loading path
-----------------------

The loading path is defined in a text file (typically ``path.txt``) located in the ``data`` folder. The file structure is as follows:

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

For tabular loading, the evolution is read from an external file:

.. code-block:: none

    #Mode
    3
    #File
    tabular_file.txt
    #Dn_init 1.
    #Dn_mini 0.01
    #prescribed_mechanical_state
    S
    0  S
    0  0  0
    #T_is_set
    0

The **#prescribed_mechanical_state** block specifies which components are controlled:

- **S**: Stress-controlled component (read from file)
- **E**: Strain-controlled component (read from file)
- **0**: Component kept constant

**#T_is_set**: 0 if temperature is constant, T if temperature is read from file.

The tabular file structure:

.. code-block:: none

    0    0.0     10   10        
    1    0.01    20   20
    2    0.02    30   30
    3    0.03    30   30
    ...

Columns: **ninc**, **time**, followed by the controlled components in order 11, 12, 22, 13, 23, 33.

If temperature is set:

.. code-block:: none

    0    0.0     293.15  10   10        
    1    0.01    294.15  20   20
    ...

Columns: **ninc**, **time**, **T**, then mechanical components.

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