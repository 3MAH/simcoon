Abaqus Integration
==================

This guide explains how to use simcoon's built-in constitutive models within Abaqus on Linux systems.

Overview
--------

Simcoon provides ready-to-use UMAT bridge files in the ``software/`` directory that allow you to use **all simcoon constitutive models** directly in Abaqus:

- ``software/umat_singleM.cpp`` - Single mechanical model (selected by material name)
- ``software/umat_singleT.cpp`` - Single thermo-mechanical model
- ``software/umat_singleM_multi.cpp`` - Multiscale mechanical model (reads from ``material.dat``)
- ``software/umat_externalM.cpp`` - Template for adding custom external UMAT in C++
- ``software/umat_externalT.cpp`` - Template for custom external thermo-mechanical UMAT

These files export an ``extern "C" umat_()`` function that Abaqus calls directly. The bridge:

1. Converts Abaqus arrays → simcoon Armadillo format (``abaqus2smart_M()``)
2. Calls simcoon's ``select_umat_M()`` to run the appropriate model based on the 5-character material name
3. Converts results back → Abaqus arrays (``smart2abaqus_M()``)

.. note::
   
   This is the **opposite** of the external plugin system documented in :doc:`../external`. 
   Here, simcoon models run inside Abaqus. In the external plugin system, external UMATs run inside simcoon's solver.

.. contents:: Table of Contents
   :local:
   :depth: 2

Prerequisites
-------------

**Linux System Requirements**

.. code-block:: bash

    # Required packages (Ubuntu/Debian)
    sudo apt-get install build-essential cmake
    sudo apt-get install libarmadillo-dev liblapack-dev libblas-dev
    
    # Required packages (RHEL/CentOS)
    sudo yum groupinstall "Development Tools"
    sudo yum install cmake armadillo-devel lapack-devel blas-devel

**simcoon Installation**

Ensure simcoon is built and installed using the provided install script:

.. code-block:: bash

    cd simcoon
    sh Install.sh -n 4
    
    # The script will:
    # - Build simcoon with CMake/Ninja
    # - Install to $CONDA_PREFIX (your conda environment)
    # - Install the Python package
    # - Optionally run tests
    
    # To skip tests:
    sh Install.sh -t -n 4

Using umat_singleM (Recommended)
--------------------------------

This is the simplest approach: the model is selected by the **material name** in your Abaqus input file (first 5 characters).

**Step 1: Compile the UMAT**

.. code-block:: bash

    cd simcoon/software
    
    # Compile to shared library (using conda environment)
    g++ -shared -fPIC -std=c++17 -O2 -o libumat_simcoon.so umat_singleM.cpp \
        -I$CONDA_PREFIX/include \
        -L$CONDA_PREFIX/lib -lsimcoon \
        -larmadillo -llapack -lblas

**Step 2: Configure Abaqus environment**

Create or edit your ``abaqus_v6.env`` file:

.. code-block:: python

    import os
    conda_prefix = os.environ.get('CONDA_PREFIX', '/path/to/conda/env')
    
    # Add library paths
    os.environ['LD_LIBRARY_PATH'] = conda_prefix + '/lib:' + os.environ.get('LD_LIBRARY_PATH', '')

**Step 3: Run Abaqus**

.. code-block:: bash

    abaqus job=mymodel user=libumat_simcoon.so

**Material definition in Abaqus input file**

The first 5 characters of the material name select the constitutive model:

.. code-block:: none

    *MATERIAL, NAME=ELISO
    *DEPVAR
    5
    *USER MATERIAL, CONSTANTS=3
    ** E, nu, alpha
    70000., 0.3, 1.2e-5

.. code-block:: none

    *MATERIAL, NAME=EPICP
    *DEPVAR
    7
    *USER MATERIAL, CONSTANTS=5
    ** E, nu, alpha, sigma_y, H
    210000., 0.3, 1.2e-5, 200., 1000.

Using umat_singleT (Thermo-mechanical)
--------------------------------------

For coupled thermo-mechanical analysis with heat generation:

.. code-block:: bash

    g++ -shared -fPIC -std=c++17 -O2 -o libumat_simcoon_T.so umat_singleT.cpp \
        -I$SIMCOON_DIR/include \
        -L$SIMCOON_DIR/lib -lsimcoon \
        -larmadillo -llapack -lblas

The thermo-mechanical version provides:

- Mechanical tangent ``ddsdde`` (∂σ/∂ε)
- Thermal stress tangent ``ddsddt`` (∂σ/∂T)
- Heat flux derivative ``drplde`` (∂r/∂ε)
- Heat capacity ``drpldt`` (∂r/∂T)
- Heat generation rate ``rpl``

Using umat_singleM_multi (Multiscale)
-------------------------------------

For multiscale homogenization models, the material definition is read from a ``data/material.dat`` file in the working directory:

.. code-block:: bash

    g++ -shared -fPIC -std=c++17 -O2 -o libumat_simcoon_multi.so umat_singleM_multi.cpp \
        -I$SIMCOON_DIR/include \
        -L$SIMCOON_DIR/lib -lsimcoon \
        -larmadillo -llapack -lblas

Create ``data/material.dat`` in your Abaqus working directory with the material definition. See the homogenization documentation for file format details.

Using umat_externalM (Custom Model)
-----------------------------------

If you want to implement a custom constitutive model in C++:

1. Copy ``software/umat_externalM.cpp`` to your working directory
2. Implement your model in the ``external_umat()`` function
3. Compile and link:

.. code-block:: bash

    g++ -shared -fPIC -std=c++17 -O2 -o libumat_custom.so umat_externalM.cpp \
        -I$SIMCOON_DIR/include \
        -L$SIMCOON_DIR/lib -lsimcoon \
        -larmadillo -llapack -lblas

The ``external_umat()`` signature receives simcoon-format data:

.. code-block:: cpp

    void external_umat(
        const vec &Etot,      // Total strain
        const vec &DEtot,     // Strain increment
        vec &sigma,           // Stress (output)
        mat &Lt,              // Tangent modulus (output)
        const mat &DR,        // Rotation increment
        const int &nprops,    // Number of properties
        const vec &props,     // Material properties
        const int &nstatev,   // Number of state variables
        vec &statev,          // State variables (input/output)
        const double &T,      // Temperature
        const double &DT,     // Temperature increment
        const double &Time,   // Step time
        const double &DTime,  // Time increment
        double &Wm,           // Mechanical work (output)
        double &Wm_r,         // Recoverable work (output)
        double &Wm_ir,        // Irrecoverable work (output)
        double &Wm_d,         // Dissipated work (output)
        const int &ndi,       // Number of direct components
        const int &nshr,      // Number of shear components
        const bool &start,    // First increment flag
        double &tnew_dt       // Suggested time step ratio (output)
    );

How It Works
------------

Architecture
~~~~~~~~~~~~

.. code-block:: none

    ┌───────────────────────────────────────────────────────────────┐
    │                        Abaqus                                 │
    │                           │                                   │
    │                           ▼                                   │
    │              calls extern "C" umat_()                         │
    └───────────────────────────────────────────────────────────────┘
                                │
                                ▼
    ┌───────────────────────────────────────────────────────────────┐
    │              umat_singleM.cpp (bridge)                        │
    │                                                               │
    │  ┌─────────────────────────────────────────────────────────┐  │
    │  │  abaqus2smart_M()                                       │  │
    │  │  - Converts stress, ddsdde, stran, dstran → Armadillo   │  │
    │  │  - Extracts temperature, properties, state variables    │  │
    │  │  - Handles 2D/3D dimensionality (ndi, nshr)             │  │
    │  └─────────────────────────────────────────────────────────┘  │
    │                           │                                   │
    │                           ▼                                   │
    │  ┌─────────────────────────────────────────────────────────┐  │
    │  │  select_umat_M(rve, ...)                                │  │
    │  │  - Reads material name (first 5 chars of cmname)        │  │
    │  │  - Dispatches to appropriate model: ELISO, EPICP, etc.  │  │
    │  └─────────────────────────────────────────────────────────┘  │
    │                           │                                   │
    │                           ▼                                   │
    │  ┌─────────────────────────────────────────────────────────┐  │
    │  │  smart2abaqus_M()                                       │  │
    │  │  - Converts sigma, Lt → stress, ddsdde arrays           │  │
    │  │  - Packs state variables + work quantities              │  │
    │  └─────────────────────────────────────────────────────────┘  │
    └───────────────────────────────────────────────────────────────┘
                                │
                                ▼
    ┌───────────────────────────────────────────────────────────────┐
    │                        Abaqus                                 │
    │              (continues with updated stress/tangent)          │
    └───────────────────────────────────────────────────────────────┘

Why extern "C" works
~~~~~~~~~~~~~~~~~~~~

The ``extern "C"`` declaration:

1. **Disables C++ name mangling** - The function is exported as ``umat_`` exactly
2. **Uses C calling convention** - Compatible with Fortran calling convention
3. **Standard types only** - ``double*``, ``int*``, ``char*`` work universally

Abaqus looks for the symbol ``umat_`` in the shared library - it doesn't care if the library was compiled from Fortran or C++.

State Variables Layout
~~~~~~~~~~~~~~~~~~~~~~

simcoon uses a specific layout for state variables (``statev`` array):

.. code-block:: none

    statev[0:nstatev_smart]  - Model-specific state variables
    statev[nstatev-4]        - Wm (total mechanical work)
    statev[nstatev-3]        - Wm_r (recoverable work)
    statev[nstatev-2]        - Wm_ir (irrecoverable work)
    statev[nstatev-1]        - Wm_d (dissipated work)

Set ``*DEPVAR`` in your input file to ``nstatev_smart + 4``.

Available Models
----------------

The following constitutive models are available through ``select_umat_M()``:

.. list-table:: Mechanical Models (umat_singleM)
   :header-rows: 1
   :widths: 12 25 30 33

   * - Code
     - Model
     - Properties (in order)
     - State Variables
   * - ELISO
     - Isotropic elasticity
     - E, ν, α
     - 1
   * - ELIST
     - Transversely isotropic elasticity
     - axis, EL, ET, νTL, νTT, GLT, αL, αT
     - 1
   * - ELORT
     - Orthotropic elasticity
     - E₁, E₂, E₃, ν₁₂, ν₁₃, ν₂₃, G₁₂, G₁₃, G₂₃, α₁, α₂, α₃
     - 1
   * - EPICP
     - Von Mises plasticity, power-law isotropic hardening
     - E, ν, α, σ_Y, k, m
     - 8 (T_init, p, EP)
   * - EPKCP
     - Von Mises, power-law isotropic + Prager kinematic
     - E, ν, α, σ_Y, k, m, kX
     - 14 (T_init, p, EP, a)
   * - EPCHA
     - Von Mises + Voce + 2× Armstrong-Frederick
     - E, ν, α, σ_Y, Q, b, C₁, D₁, C₂, D₂
     - 33
   * - EPHIL / EPTRI
     - Hill yield + power-law isotropic hardening
     - E, ν, α, σ_Y, k, m, F, G, H, L, M, N
     - 8 (T_init, p, EP)
   * - EPHAC
     - Cubic elasticity + Hill + Voce + 2× AF
     - E, ν, G, α, σ_Y, Q, b, C₁, D₁, C₂, D₂, F, G, H, L, M, N
     - 33
   * - EPANI
     - Cubic elasticity + anisotropic yield + Voce + 2× AF
     - E, ν, G, α, σ_Y, Q, b, C₁, D₁, C₂, D₂, P₁₁..P₆₆ (9)
     - 33
   * - EPDFA
     - Cubic elasticity + DFA yield + Voce + 2× AF
     - E, ν, G, α, σ_Y, Q, b, C₁, D₁, C₂, D₂, F, G, H, L, M, N, K
     - 33
   * - EPCHG
     - Generic Chaboche (selectable yield, N iso/kin terms)
     - E, ν, G, α, σ_Y, N_iso, N_kin, criteria, (Q,b)×N, (C,D)×N, crit. params
     - 33
   * - EPHIN
     - N Hill yield surfaces
     - E, ν, α, N, per surface: σ_Y, k, m, F, G, H, L, M, N
     - 1 + 7N
   * - SMAUT
     - SMA unified model
     - See SMA documentation
     - 24
   * - SMANI
     - SMA anisotropic model
     - See SMA documentation
     - 24
   * - LLDM0
     - Lemaitre-Ladeveze-Dufailly damage
     - See header documentation
     - 9
   * - ZENER
     - Kelvin viscoelastic (single branch)
     - E₀, ν₀, α, E₁, ν₁, ηB₁, ηS₁
     - 14
   * - ZENNK
     - Kelvin viscoelastic (N branches)
     - E₀, ν₀, α, N, per branch: Eᵢ, νᵢ, ηBᵢ, ηSᵢ
     - 7 + 7N
   * - PRONK
     - Prony series viscoelastic (generalized Maxwell)
     - E₀, ν₀, α, N, per branch: Eᵢ, νᵢ, ηBᵢ, ηSᵢ
     - 7 + 7N
   * - MODUL
     - Composable modular UMAT
     - self-describing stream (see :mod:`simcoon.modular`)
     - model-dependent

Several of these names are served by the modular engine through
props-translating adapters — identical usage and results; see
:doc:`umat_catalog` for per-name status and state-variable layout notes.

.. list-table:: Micromechanics Models
   :header-rows: 1
   :widths: 12 40 48

   * - Code
     - Model
     - Description
   * - MIHEN
     - Mori-Tanaka (Eshelby)
     - Mean-field homogenization with Eshelby tensor
   * - MIMTN
     - Mori-Tanaka with N phases
     - Multi-phase Mori-Tanaka
   * - MISCN
     - Self-consistent (N phases)
     - Self-consistent homogenization
   * - MIPLN
     - Periodic layered
     - Layered composite homogenization

For micromechanics models, use ``umat_singleM_multi.cpp`` with a ``data/material.dat`` file.

Troubleshooting
---------------

**Compilation errors**

.. code-block:: bash

    # Check include paths
    ls $SIMCOON_DIR/include/simcoon/
    
    # Check library exists
    ls $SIMCOON_DIR/lib/libsimcoon.*

**Runtime library errors**

.. code-block:: bash

    # Check library path
    echo $LD_LIBRARY_PATH
    
    # Check dependencies of your UMAT
    ldd libumat_simcoon.so
    
    # Look for missing libraries
    ldd libumat_simcoon.so | grep "not found"

**Symbol not found errors**

.. code-block:: bash

    # Verify umat_ symbol is exported
    nm -D libumat_simcoon.so | grep umat_
    # Should show: T umat_

**Convergence issues**

- Check that ``*DEPVAR`` matches ``nstatev_smart + 4``
- Verify material properties are in correct units
- Start with small load increments
- Check the tangent modulus is properly computed (symmetric, positive definite for elastic)

**Debugging**

Add debug output in the bridge code:

.. code-block:: cpp

    #include <fstream>
    
    // In umat_() function, after abaqus2smart_M():
    std::ofstream debug("umat_debug.log", std::ios::app);
    debug << "=== Increment " << kinc << " ===" << std::endl;
    debug << "Material: " << umat_name << std::endl;
    debug << "T = " << rve_sv_M->T << ", DT = " << rve_sv_M->DT << std::endl;
    debug << "Etot: " << rve_sv_M->Etot.t();
    debug << "DEtot: " << rve_sv_M->DEtot.t();
    debug.close();

Complete Example
----------------

Here's a complete example for a tensile test on a plasticity model:

.. code-block:: none

    *HEADING
    Tensile test with simcoon EPICP model
    **
    *NODE
    1, 0.0, 0.0, 0.0
    2, 1.0, 0.0, 0.0
    3, 1.0, 1.0, 0.0
    4, 0.0, 1.0, 0.0
    5, 0.0, 0.0, 1.0
    6, 1.0, 0.0, 1.0
    7, 1.0, 1.0, 1.0
    8, 0.0, 1.0, 1.0
    **
    *ELEMENT, TYPE=C3D8, ELSET=ALL
    1, 1, 2, 3, 4, 5, 6, 7, 8
    **
    *SOLID SECTION, ELSET=ALL, MATERIAL=EPICP
    **
    *MATERIAL, NAME=EPICP
    *DEPVAR
    12
    *USER MATERIAL, CONSTANTS=5
    ** E, nu, alpha, sigma_y, H
    210000., 0.3, 0., 200., 1000.
    **
    *BOUNDARY
    1, 1, 3
    2, 2, 3
    3, 3
    4, 1
    4, 3
    5, 1, 2
    6, 2
    8, 1
    **
    *STEP, NLGEOM=NO
    *STATIC
    0.1, 1.0, 1e-5, 0.1
    **
    *BOUNDARY
    5, 3, 3, 0.01
    6, 3, 3, 0.01
    7, 3, 3, 0.01
    8, 3, 3, 0.01
    **
    *OUTPUT, FIELD
    *NODE OUTPUT
    U, RF
    *ELEMENT OUTPUT
    S, E, SDV
    **
    *END STEP
