Ansys Integration
=================

This guide explains how to use simcoon's built-in constitutive models within Ansys MAPDL on Linux systems.

.. note::

   Ansys USERMAT integration is more complex than Abaqus due to differences in the interface design.
   Simcoon provides ``software/usermat_singleM.cpp`` as the Ansys equivalent of ``umat_singleM.cpp``.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
--------

The Ansys ``USERMAT`` subroutine interface differs from Abaqus in several ways:

- **Voigt notation**: Ansys uses (11, 22, 33, 12, 23, 13) vs simcoon's (11, 22, 33, 12, 13, 23)
- **Finite strain handling**: Built-in Jaumann framework for rate-form behaviours
- **No time step control**: Cannot request smaller time steps from within USERMAT
- **State variable initialization**: Cannot set non-zero initial values
- **Model selection**: Uses numeric model code in props[0] instead of material name

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

Using usermat_singleM.cpp
-------------------------

The ``software/usermat_singleM.cpp`` file provides the same functionality as the Abaqus bridge, but adapted for Ansys.

**Step 1: Compile the USERMAT**

.. code-block:: bash

    cd simcoon/software
    
    # Compile to shared library (using conda environment)
    g++ -shared -fPIC -std=c++17 -O2 -o libusermat_simcoon.so usermat_singleM.cpp \
        -I$CONDA_PREFIX/include \
        -L$CONDA_PREFIX/lib -lsimcoon \
        -larmadillo -llapack -lblas

**Step 2: Set up library path**

Ansys searches for libraries in directories specified by ``ANS_USER_PATH``:

.. code-block:: bash

    export ANS_USER_PATH=/path/to/your/usermat/directory
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$ANS_USER_PATH:$LD_LIBRARY_PATH

**Step 3: Define material in MAPDL**

Unlike Abaqus where the model is selected by material name, Ansys uses a **model code** as the first material property:

.. code-block:: none

    ! EPICP model example (model code 5)
    ! Properties: code, E, nu, alpha, sigma_y, H
    TB,USER,1,1,6
    TBDATA,1,5         ! Model code for EPICP
    TBDATA,2,210000    ! E
    TBDATA,3,0.3       ! nu
    TBDATA,4,0.0       ! alpha
    TBDATA,5,200       ! sigma_y
    TBDATA,6,1000      ! H
    
    ! Define state variables (nstatev_smart + 4)
    TB,STATE,1,,12

Model Codes
-----------

.. list-table:: Ansys Model Codes
   :header-rows: 1
   :widths: 10 15 35 40

   * - Code
     - Name
     - Model
     - Additional Properties
   * - 2
     - ELISO
     - Isotropic elasticity
     - E, ν, α
   * - 3
     - ELIST
     - Transversely isotropic elasticity
     - E₁, E₂, ν₁₂, ν₂₃, G₁₂, α₁, α₂
   * - 4
     - ELORT
     - Orthotropic elasticity
     - E₁, E₂, E₃, ν₁₂, ν₁₃, ν₂₃, G₁₂, G₁₃, G₂₃, α₁, α₂, α₃
   * - 5
     - EPICP
     - Isotropic plasticity (isotropic hardening)
     - E, ν, α, σ_y, H
   * - 6
     - EPKCP
     - Kinematic + isotropic hardening
     - E, ν, α, σ_y, H, C, γ
   * - 7
     - EPCHA
     - Chaboche cyclic plasticity
     - E, ν, α, σ_y, Q, b, C₁, γ₁, ...
   * - 8
     - SMAUT
     - SMA unified model
     - See SMA documentation
   * - 10
     - LLDM0
     - Lemaitre-Chaboche damage
     - E, ν, α, σ_y, H, S, s, D_c
   * - 11
     - ZENER
     - Zener viscoelastic (single branch)
     - E₀, E₁, η
   * - 12
     - ZENNK
     - Zener viscoelastic (N branches)
     - E₀, E₁, η₁, E₂, η₂, ...
   * - 13
     - PRONK
     - Prony series viscoelastic
     - G₀, K₀, g₁, τ₁, k₁, τ'₁, ...
   * - 17
     - EPHIL
     - Hill anisotropic plasticity
     - E, ν, α, σ_y, H, F, G, H, L, M, N
   * - 18
     - EPHAC
     - Hill + Chaboche
     - E, ν, α, σ_y, Q, b, C₁, γ₁, F, G, H, L, M, N

How It Works
------------

Architecture
~~~~~~~~~~~~

.. code-block:: none

    ┌───────────────────────────────────────────────────────────────┐
    │                     Ansys MAPDL                               │
    │                           │                                   │
    │                           ▼                                   │
    │              calls extern "C" usermat_()                      │
    └───────────────────────────────────────────────────────────────┘
                                │
                                ▼
    ┌───────────────────────────────────────────────────────────────┐
    │              usermat_singleM.cpp (bridge)                     │
    │                                                               │
    │  ┌─────────────────────────────────────────────────────────┐  │
    │  │  ansys2smart_M()                                        │  │
    │  │  - Converts stress, dsdePl, strain, dstrain → Armadillo │  │
    │  │  - Remaps Voigt indices (swaps 4 ↔ 5 for 3D)            │  │
    │  │  - Extracts model code from props[0]                    │  │
    │  └─────────────────────────────────────────────────────────┘  │
    │                           │                                   │
    │                           ▼                                   │
    │  ┌─────────────────────────────────────────────────────────┐  │
    │  │  select_umat_M(rve, ...)                                │  │
    │  │  - Looks up model name from code (e.g., 5 → "EPICP")    │  │
    │  │  - Dispatches to appropriate model                      │  │
    │  └─────────────────────────────────────────────────────────┘  │
    │                           │                                   │
    │                           ▼                                   │
    │  ┌─────────────────────────────────────────────────────────┐  │
    │  │  smart2ansys_M()                                        │  │
    │  │  - Converts sigma, Lt → stress, dsdePl arrays           │  │
    │  │  - Remaps Voigt indices back to Ansys convention        │  │
    │  └─────────────────────────────────────────────────────────┘  │
    └───────────────────────────────────────────────────────────────┘
                                │
                                ▼
    ┌───────────────────────────────────────────────────────────────┐
    │                     Ansys MAPDL                               │
    │              (continues with updated stress/tangent)          │
    └───────────────────────────────────────────────────────────────┘

Voigt Notation
~~~~~~~~~~~~~~

simcoon and Abaqus use the same Voigt notation, but Ansys differs in the shear component ordering:

.. list-table:: Voigt Index Comparison
   :header-rows: 1
   :widths: 15 25 30 30

   * - Index
     - Component
     - simcoon / Abaqus
     - Ansys
   * - 0
     - Normal 11
     - σ₁₁, ε₁₁
     - σ₁₁, ε₁₁
   * - 1
     - Normal 22
     - σ₂₂, ε₂₂
     - σ₂₂, ε₂₂
   * - 2
     - Normal 33
     - σ₃₃, ε₃₃
     - σ₃₃, ε₃₃
   * - 3
     - Shear 12
     - σ₁₂, γ₁₂
     - σ₁₂, γ₁₂
   * - 4
     - Shear 13/23
     - σ₁₃, γ₁₃
     - σ₂₃, γ₂₃
   * - 5
     - Shear 23/13
     - σ₂₃, γ₂₃
     - σ₁₃, γ₁₃

The ``usermat_singleM.cpp`` bridge automatically swaps indices 4 and 5 during conversion.

State Variables Layout
~~~~~~~~~~~~~~~~~~~~~~

Same as Abaqus - simcoon uses a specific layout for state variables:

.. code-block:: none

    ustatev[0:3]              - Work quantities (Wm, Wm_r, Wm_ir, Wm_d)
    ustatev[4:nStatev-1]      - Model-specific state variables

Set the number of state variables in MAPDL with ``TB,STATE,mat_id,,nstatev``.

Limitations and Workarounds
---------------------------

Based on MFront's experience with Ansys:

1. **No external state variables**: Temperature must be passed through material properties if needed beyond the standard interface.

2. **No non-zero initial state variables**: If your model requires non-zero initial values, consider:
   
   - Initializing in the first increment based on material properties
   - Using a separate initialization step

3. **No time step control**: Unlike Abaqus where ``pnewdt`` can request smaller steps, Ansys USERMAT cannot control the time step. Use sufficiently small increments for nonlinear problems.

4. **Modelling hypothesis ambiguity**: Plane strain and axisymmetric cases have identical tensor dimensions. If your model treats them differently, consider using a material property flag to distinguish them.

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
    echo $ANS_USER_PATH
    
    # Check dependencies of your USERMAT
    ldd libusermat_simcoon.so
    
    # Look for missing libraries
    ldd libusermat_simcoon.so | grep "not found"

**Symbol not found errors**

.. code-block:: bash

    # Verify usermat_ symbol is exported
    nm -D libusermat_simcoon.so | grep usermat_
    # Should show: T usermat_

**Convergence issues**

- Check that ``TB,STATE`` count matches ``nstatev_smart + 4``
- Verify material properties are in correct units
- Start with small load increments
- The first property must be the model code

References
----------

- `MFront Ansys Interface Documentation <https://thelfer.github.io/tfel/web/ansys.html>`_
- Ansys MAPDL User Material Subroutine Guide
- Belytschko, T. *Nonlinear Finite Elements for Continua and Structures*. Wiley, 2000.
