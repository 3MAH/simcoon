/**
 * @file simcoon_doxygen_groups.hpp
 * @brief Doxygen group definitions for Simcoon library documentation
 * 
 * This file defines all the Doxygen groups used to organize the API documentation.
 * It should not be included in code - it is only for documentation purposes.
 */

#pragma once

//=============================================================================
// Continuum Mechanics - Functions
//=============================================================================

/**
 * @defgroup constitutive Constitutive Functions
 * @brief Utility functions for constitutive modeling
 * 
 * Functions for building stiffness and compliance tensors, identity tensors
 * in Voigt notation, and elastic predictions.
 */

/**
 * @defgroup contimech Continuum Mechanics Utilities
 * @brief General continuum mechanics tensor operations
 * 
 * Functions for tensor transformations, Voigt notation conversions,
 * and common continuum mechanics operations.
 */

/**
 * @defgroup kinematics Kinematics Functions
 * @brief Finite strain kinematics functions
 * 
 * Functions for computing deformation gradients, stretch tensors,
 * rotation tensors, and strain measures.
 */

/**
 * @defgroup hyperelastic Hyperelastic Functions
 * @brief Hyperelastic material model functions
 * 
 * Functions for strain energy derivatives, stress computations,
 * and tangent moduli for hyperelastic models.
 */

/**
 * @defgroup stress Stress Functions
 * @brief Stress transformation functions
 * 
 * Functions for stress transformations between different configurations
 * and stress measures.
 */

/**
 * @defgroup criteria Yield Criteria Functions
 * @brief Yield and failure criteria functions
 * 
 * Functions for computing yield functions and their derivatives
 * for various yield criteria (von Mises, Hill, etc.).
 */

/**
 * @defgroup damage Damage Functions
 * @brief Damage mechanics functions
 * 
 * Functions for damage evolution and degradation modeling.
 */

/**
 * @defgroup objective_rates Objective Rates Functions
 * @brief Objective stress rate functions
 * 
 * Functions for computing objective stress rates
 * (Jaumann, Green-Naghdi, Truesdell, etc.).
 */

/**
 * @defgroup recovery_props Recovery Properties Functions
 * @brief Material property recovery functions
 * 
 * Functions for recovering material properties from experimental data
 * or stiffness tensors.
 */

/**
 * @defgroup transfer Transfer Functions
 * @brief Configuration transfer functions
 * 
 * Functions for transferring quantities between reference
 * and current configurations.
 */

/**
 * @defgroup derivatives Tensor Derivatives
 * @brief Tensor derivative functions
 * 
 * Functions for computing derivatives of tensor invariants,
 * determinants, and inverses.
 */

//=============================================================================
// Continuum Mechanics - Homogenization
//=============================================================================

/**
 * @defgroup eshelby Eshelby Tensor
 * @brief Eshelby tensor computation functions
 * 
 * Functions for computing Eshelby tensors for various
 * inclusion geometries (ellipsoids, cylinders, layers).
 */

/**
 * @defgroup homogenization Homogenization Functions
 * @brief Mean-field homogenization functions
 * 
 * Functions for computing effective properties using
 * Mori-Tanaka, self-consistent, and other schemes.
 */

//=============================================================================
// Continuum Mechanics - Material
//=============================================================================

/**
 * @defgroup material Material Properties
 * @brief Material characterization functions
 * 
 * Functions and classes for material property management,
 * ODF/PDF processing, and crystallography.
 */

//=============================================================================
// Continuum Mechanics - Micromechanics
//=============================================================================

/**
 * @defgroup micromechanics Micromechanics Functions
 * @brief Micromechanical modeling functions
 * 
 * Functions for multiphase material modeling and
 * micromechanical schemes.
 */

//=============================================================================
// UMAT - Constitutive Models
//=============================================================================

/**
 * @defgroup umat_finite Finite Strain UMAT
 * @brief Finite strain constitutive models
 * 
 * User material subroutines for hyperelastic models
 * (Neo-Hookean, Mooney-Rivlin, etc.).
 */

/**
 * @defgroup umat_mechanical Mechanical UMAT
 * @brief Small strain mechanical constitutive models
 * 
 * User material subroutines for elasticity, plasticity,
 * viscoelasticity, and damage models.
 */

/**
 * @defgroup umat_thermomechanical Thermomechanical UMAT
 * @brief Coupled thermomechanical constitutive models
 * 
 * User material subroutines for temperature-dependent
 * material behavior.
 */

//=============================================================================
// Simulation
//=============================================================================

/**
 * @defgroup solver Solver Functions
 * @brief Simulation solver functions
 * 
 * Functions for solving thermomechanical boundary value problems.
 */

/**
 * @defgroup identification Identification Functions
 * @brief Parameter identification functions
 * 
 * Functions for material parameter identification using
 * optimization algorithms.
 */

/**
 * @defgroup maths Mathematical Functions
 * @brief Mathematical utility functions
 * 
 * Functions for rotations, random number generation,
 * statistics, and numerical solving.
 */

/**
 * @defgroup phase Phase Management
 * @brief Phase and state variable management
 * 
 * Classes and functions for managing material phases
 * and internal state variables.
 */

/**
 * @defgroup geometry Geometry Functions
 * @brief Inclusion geometry functions
 * 
 * Classes for defining inclusion geometries
 * (ellipsoids, cylinders, layers).
 */
