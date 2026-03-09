======================
Mathematical Functions
======================

Mathematical utility functions for rotations, random number generation,
statistics, and numerical solving.

.. toctree::
   :maxdepth: 2

   rotation

Overview
========

The Maths module provides essential mathematical utilities for continuum mechanics
simulations:

- **Rotation**: Comprehensive 3D rotation tools including the new ``Rotation`` class
  and Voigt notation transformations (see :doc:`rotation`)
- **Random**: Random number generation utilities
- **Statistics**: Statistical functions
- **Numerical Solvers**: Root finding and optimization

Random Number Generation
========================

.. doxygengroup:: random
   :project: simcoon
   :content-only:

Statistics
==========

.. doxygengroup:: stats
   :project: simcoon
   :content-only:

Numerical Solvers
=================

.. doxygengroup:: solvers
   :project: simcoon
   :content-only:
