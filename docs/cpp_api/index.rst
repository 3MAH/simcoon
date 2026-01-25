====================
C++ API Reference
====================

This section contains the C++ API documentation for Simcoon, generated from the source code using Doxygen and integrated via Breathe.

Simcoon's C++ library provides high-performance implementations for:

- **Continuum Mechanics**: Tensor operations, strain/stress transformations, kinematics
- **Constitutive Models**: Hyperelastic, viscoelastic, and viscoplastic material models
- **Homogenization**: Mean-field homogenization methods for composite materials
- **Simulation Tools**: Solvers, identification algorithms, and utilities

.. note::
   The C++ API is designed for advanced users who need direct access to the 
   computational core. For most use cases, the Python bindings provide a more 
   convenient interface.

Namespace Overview
==================

All Simcoon C++ functions and classes are contained within the ``simcoon`` namespace:

.. code-block:: cpp

   #include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>
   
   using namespace simcoon;
   
   // Example: compute isochoric invariants
   arma::mat F = arma::randu(3,3);
   arma::mat b = L_Cauchy_Green(F);
   double J = arma::det(F);
   arma::vec I_bar = isochoric_invariants(b, J);

.. toctree::
   :maxdepth: 2
   :caption: Modules:

   continuum_mechanics/index
   simulation/index
   homogenization/index
   umat/index
