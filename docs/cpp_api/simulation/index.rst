==========
Simulation
==========

The Simulation module provides tools for running simulations, including solvers,
identification algorithms, and mathematical utilities.

.. toctree::
   :maxdepth: 2
   :caption: Submodules:

   maths
   rotation
   solver
   identification
   phase

Overview
========

This module contains simulation and numerical tools:

- **Maths**: Mathematical utilities (random numbers, statistics, solvers)
- **Rotation**: Comprehensive 3D rotation tools with ``Rotation`` class and Voigt notation support
- **Solver**: Material point simulation solvers
- **Identification**: Parameter identification algorithms
- **Phase**: Phase management and properties

What's New in Simcoon 2.0
=========================

The **Rotation module** has been significantly enhanced with:

- New ``Rotation`` class inspired by ``scipy.spatial.transform.Rotation``
- Unit quaternion internal representation for numerical stability
- Support for multiple Euler angle conventions (zxz, zyz, xyz, etc.)
- Intrinsic and extrinsic rotation modes
- SLERP interpolation for smooth orientation transitions
- Direct Voigt notation operations (stress, strain, stiffness, compliance)
- Full Python bindings with NumPy integration

See :doc:`rotation` for complete documentation.
