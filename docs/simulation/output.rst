Configure the output
================================

The output options for the solver are configured through an ``output.dat`` file located in the ``data/`` directory. This file allows you to customize what quantities are written to the result file, including the type of strain and stress measures, rotation information, tangent modulus, temperature, and internal state variables.

File structure
--------------

The ``output.dat`` file should be placed in the ``data/`` subdirectory relative to your simulation working directory. The complete structure of the file is as follows:

.. code-block:: none

    #Output_values
    strain_type <value>
    nb_strain   <number>
    <strain component indices...>
    stress_type <value>
    nb_stress   <number>
    <stress component indices...>

    Rotation_type <value>
    Tangent_type  <value>
    T             <value>

    Number_of_wanted_internal_variables <n|all>
    [#number from <start> to <end>]

    #Block #type_1_N_2_T    #every
    <block_num> <output_type> <frequency>
    ...

Basic example
-------------

Here is a basic example of an ``output.dat`` file:

.. code-block:: none

    #Output_values
    strain_type 0
    nb_strain   6
    0   1   2   3   4   5
    stress_type	4
    nb_stress   6
    0   1   2   3   4   5

    Rotation_type	0
    Tangent_type	0
    T   1

    Number_of_wanted_internal_variables	0

    #Block #type_1_N_2_T    #every
    1      1                1

This configuration outputs all 6 components of Green-Lagrange strain and Cauchy stress, along with temperature, at every increment of block 1.

.. note::

   In the case of small deformations (Control_type = 1, no NLGEOM), all strain and stress measures reduce to the same quantities. The Green-Lagrange strain becomes equivalent to the infinitesimal strain tensor, and all stress measures become equivalent to the Cauchy stress. Therefore, strain_type = 0 (Green-Lagrange) and stress_type = 4 (Cauchy) are used by default for small deformation problems.

Strain output options
---------------------

**strain_type** controls which strain measure is output. The available options are:

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Value
     - Strain Type
     - Description
   * - 0
     - Green-Lagrange
     - :math:`\mathbf{E} = \frac{1}{2}(\mathbf{F}^T\mathbf{F} - \mathbf{I})` in Voigt notation (6 components)
   * - 1
     - Biot
     - :math:`\sqrt{2\mathbf{E}+\mathbf{I}} - \mathbf{I}` (6 components)
   * - 2
     - Deformation Gradient
     - Full :math:`\mathbf{F}` matrix (9 components) as row-wise vector
   * - 3
     - Logarithmic (Eulerian)
     - Small strain tensor :math:`\boldsymbol{\varepsilon}` (6 components)
   * - 4
     - Isochoric Principal Stretches
     - :math:`\bar{\lambda}_i` (3 principal stretches from :math:`\mathbf{V}` decomposition)

**nb_strain** specifies the number of strain components to output.

The following line contains space-separated component indices. For Voigt notation (types 0, 1, 3), the indices are:

- 0: :math:`\varepsilon_{11}` (or :math:`E_{11}`)
- 1: :math:`\varepsilon_{22}`
- 2: :math:`\varepsilon_{33}`
- 3: :math:`\varepsilon_{12}` (or :math:`\gamma_{12}`)
- 4: :math:`\varepsilon_{13}`
- 5: :math:`\varepsilon_{23}`

For full tensor notation (type 2), the indices are 0-8 for the 9 components in row-wise order.

For isochoric principal stretches (type 4), the indices are 0-2 for the 3 principal stretches.

Stress output options
---------------------

**stress_type** controls which stress measure is output. The available options are:

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Value
     - Stress Type
     - Description
   * - 0
     - 2nd Piola-Kirchhoff (PKII)
     - :math:`\mathbf{S}` - Reference configuration (6 components)
   * - 1
     - Nominal Stress
     - :math:`\mathbf{P}^T` - Transpose of 1st Piola-Kirchhoff (9 components)
   * - 2
     - 1st Piola-Kirchhoff (PKI)
     - :math:`\mathbf{P}` - Unsymmetric (9 components)
   * - 3
     - Kirchhoff
     - :math:`\boldsymbol{\tau} = J\boldsymbol{\sigma}` (6 components)
   * - 4
     - Cauchy
     - :math:`\boldsymbol{\sigma}` - True stress (6 components)

**nb_stress** specifies the number of stress components to output.

The following line contains space-separated component indices. For Voigt notation (types 0, 3, 4), the indices are 0-5. For full tensor notation (types 1, 2), the indices are 0-8.

Rotation output options
-----------------------

**Rotation_type** controls what rotation information is output:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Value
     - Description
   * - 0
     - No rotation output
   * - 1
     - Rotation matrix :math:`\mathbf{R}` (9 components) + corotational basis vectors :math:`\mathbf{g}_i` (9 components)
   * - 2
     - Rotation increment :math:`\Delta\mathbf{R}` (9 components)
   * - 3
     - Rotation matrix :math:`\mathbf{R}` (9 components) + rotation increment :math:`\Delta\mathbf{R}` (9 components)

Tangent modulus output
----------------------

**Tangent_type** controls whether the tangent modulus is output:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Value
     - Description
   * - 0
     - No tangent modulus output
   * - 1
     - Output tangent modulus :math:`\mathbf{L}_t` (36 components for mechanical, additional components for thermomechanical)

For thermomechanical problems (Loading_type = 2), additional coupling terms are output: :math:`\partial\mathbf{S}/\partial T` (6 components), :math:`\partial r/\partial\mathbf{E}` (6 components), and :math:`\partial r/\partial T` (1 component).

Temperature output
------------------

**T** controls whether temperature information is output:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Value
     - Description
   * - 0
     - No temperature output
   * - 1
     - Output temperature :math:`T`, heat flux :math:`Q`, and heat source :math:`r`

Internal state variables
------------------------

**Number_of_wanted_internal_variables** controls which internal state variables are output. There are three options:

1. **0**: No internal state variables are output
2. **all** (or **All** or **ALL**): Output all internal state variables
3. **n** (a positive integer): Output n specific selections of state variables

If you specify a positive integer n, you need to provide n selection lines. Each selection can be either:

- A single state variable: ``#<num> <index>``
- A range of state variables: ``#<num> from <start> to <end>``

Example with specific state variables:

.. code-block:: none

    Number_of_wanted_internal_variables	3
    #1 5
    #2 from 10 to 15
    #3 20

This outputs state variable 5, state variables 10 through 15, and state variable 20.

Block output control
--------------------

The last section controls the output frequency for each simulation block:

.. code-block:: none

    #Block #type_1_N_2_T    #every
    <block_num> <type> <frequency>

For each block, you specify:

- **block_num**: The block number
- **type**: The output type (1 for increment-based, 2 for time-based)
- **frequency**: The output frequency

.. list-table::
   :header-rows: 1
   :widths: 10 20 70

   * - Type
     - Meaning
     - Frequency
   * - 1
     - Every N increments
     - Integer N = output every N increments
   * - 2
     - Time-based
     - Float = output at time intervals

Example for multiple blocks:

.. code-block:: none

    #Block #type_1_N_2_T    #every
    1      1                1
    2      1                5
    3      2                0.1

This outputs every increment for block 1, every 5 increments for block 2, and at time intervals of 0.1 for block 3.

Hyperelasticity example
-----------------------

Here is an example for hyperelastic simulations where you want to output isochoric principal stretches and nominal stress:

.. code-block:: none

    #Output_values
    strain_type 4
    nb_strain   3
    0   1   2
    stress_type	1
    nb_stress   9
    0   1   2   3   4   5   6   7   8

    Rotation_type	0
    Tangent_type	0
    T   1

    Number_of_wanted_internal_variables	0

    #Block #type_1_N_2_T    #every
    1      1                1

This configuration outputs the 3 isochoric principal stretches and all 9 components of the nominal stress tensor.

Default values
--------------

If the ``output.dat`` file is not present, the solver uses the following default values:

- **strain_type**: 0 (Green-Lagrange)
- **stress_type**: 4 (Cauchy)
- **nb_strain**: 6 (all components)
- **nb_stress**: 6 (all components)
- **Rotation_type**: 0 (no rotation output)
- **Tangent_type**: 0 (no tangent modulus)
- **T**: 1 (output temperature)
- **Number_of_wanted_internal_variables**: 0 (no state variables)
- **type**: 1 (increment-based for all blocks)
- **frequency**: 1 (every increment)

Result file format
------------------

The result file contains columns in the following order:

1. Block number
2. Cycle number
3. Step number
4. Increment number
5. Time
6. Temperature, heat flux Q, heat source r (if T = 1)
7. Strain components (based on strain_type and indices)
8. Stress components (based on stress_type and indices)
9. Rotation data (if Rotation_type > 0)
10. Tangent modulus (if Tangent_type = 1)
11. Mechanical energy components: :math:`W_m` (4 values)
12. Thermal energy components: :math:`W_t` (3 values, only for thermomechanical problems)
13. Selected internal state variables (if Number_of_wanted_internal_variables > 0 or all)