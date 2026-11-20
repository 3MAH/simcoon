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

- **Rotation**: the ``Rotation`` class and the Voigt rotation helpers (see :doc:`rotation`)
- **Random numbers and statistics**: random draws, probability distributions, peak functions
- **Local solvers**: Newton and Fischer-Burmeister solvers of the return-mapping systems,
  Lagrange penalty functions, quadratic roots
- **Fixed-size kernels**: 3x3 spectral decompositions and fixed-order products used in
  the per-point loops

.. _cpp-maths-functions:

Functions
=========

.. doxygengroup:: maths
   :project: simcoon
   :content-only:
