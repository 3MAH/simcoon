FEA Integration
===============

This section explains how to use simcoon's built-in constitutive models within commercial Finite Element Analysis (FEA) software on Linux systems.

.. note::
   
   This is the **opposite** of the external plugin system documented in :doc:`../external`. 
   Here, simcoon models run inside commercial FEA software. In the external plugin system, external UMATs run inside simcoon's solver.

Overview
--------

Simcoon provides ready-to-use bridge files in the ``software/`` directory that allow you to use **all simcoon constitutive models** directly in commercial FEA software:

**For Abaqus:**

- ``software/umat_singleM.cpp`` - Single mechanical model (selected by material name)
- ``software/umat_singleT.cpp`` - Single thermo-mechanical model
- ``software/umat_singleM_multi.cpp`` - Multiscale mechanical model
- ``software/umat_externalM.cpp`` - Template for custom external UMAT

**For Ansys:**

- ``software/usermat_singleM.cpp`` - Single mechanical model (selected by model code)

The bridge files:

1. Export an ``extern "C"`` function (``umat_()`` or ``usermat_()``) that the FEA software calls
2. Convert FEA arrays ↔ simcoon Armadillo format
3. Call simcoon's ``select_umat_M()`` to run the appropriate constitutive model
4. Convert results back to FEA format

Supported FEA Software
----------------------

.. toctree::
   :maxdepth: 2

   abaqus
   ansys

Key Differences
---------------

.. list-table:: Comparison of FEA Interfaces
   :header-rows: 1
   :widths: 25 35 40

   * - Feature
     - Abaqus
     - Ansys
   * - Bridge file
     - ``umat_singleM.cpp``
     - ``usermat_singleM.cpp``
   * - Entry point
     - ``umat_()``
     - ``usermat_()``
   * - Model selection
     - Material name (first 5 chars)
     - Model code in props[0]
   * - Voigt notation
     - (11,22,33,12,13,23)
     - (11,22,33,12,23,13)
   * - Time step control
     - Yes (pnewdt)
     - No
   * - Thermo-mechanical
     - ``umat_singleT.cpp``
     - Not yet implemented

Available Models
----------------

All simcoon constitutive models are available through the FEA bridges:

**Elasticity:**
ELISO (isotropic), ELIST (transversely isotropic), ELORT (orthotropic)

**Plasticity:**
EPICP (isotropic hardening), EPKCP (kinematic), EPCHA (Chaboche),
EPHIL/EPTRI (Hill), EPHAC (Hill-Chaboche), EPANI (anisotropic-Chaboche),
EPDFA (DFA-Chaboche), EPCHG (generic Chaboche), EPHIN (N Hill surfaces)

**Viscoelasticity:**
ZENER (Kelvin, single branch), ZENNK (Kelvin, N branches), PRONK (Prony series)

**Shape Memory Alloys:**
SMAUT (unified), SMANI (anisotropic)

**Damage:**
LLDM0 (Lemaitre-Ladeveze-Dufailly)

**Micromechanics:**
MIHEN (Mori-Tanaka), MIMTN (multi-phase), MISCN (self-consistent), MIPLN (layered)

**Composable:**
MODUL (modular UMAT — elasticity + any combination of plasticity,
viscoelasticity and damage mechanisms; see :mod:`simcoon.modular`)

Many legacy names are now served by the modular engine through
props-translating adapters (identical calling convention and results): see
:doc:`umat_catalog` for the complete per-name status, props layouts and
state-variable notes.

See the individual Abaqus and Ansys pages for complete property and state variable specifications.

