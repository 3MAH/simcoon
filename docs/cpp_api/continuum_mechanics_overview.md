# Continuum Mechanics Module Overview

## Introduction

The Continuum Mechanics module provides the fundamental building blocks for modeling the mechanical behavior of materials under various loading conditions. It encompasses constitutive relations, kinematic descriptions, stress transformations, homogenization methods, and material models (UMAT) for finite element analysis.

## Module Organization

The Continuum Mechanics module is organized into several submodules:

### 1. **Functions** - Core Computational Operations

Mathematical functions and tensor operations that form the foundation of continuum mechanics calculations:

- **constitutive.hpp** - Construction of elastic stiffness and compliance tensors
- **contimech.hpp** - Tensor manipulations, Voigt notation conversions, stress/strain invariants
- **kinematics.hpp** - Finite strain kinematics (deformation gradients, stretch tensors)
- **stress.hpp** - Stress measure transformations (Cauchy, Piola-Kirchhoff, Biot)
- **hyperelastic.hpp** - Hyperelastic material functions (strain energy, tangent moduli)
- **criteria.hpp** - Yield and failure criteria (von Mises, Hill, Drucker-Prager, Tresca)
- **damage.hpp** - Damage mechanics functions (Weibull, Kachanov, Miner, Manson-Coffin)
- **objective_rates.hpp** - Objective stress rates (Jaumann, Green-Naghdi, Truesdell, logarithmic)
- **recovery_props.hpp** - Material property recovery from stiffness/compliance tensors
- **transfer.hpp** - Conversions between tensor and Voigt notations
- **derivatives.hpp** - Tensor derivatives for sensitivity analysis
- **natural_basis.hpp** - Natural basis transformations for curvilinear coordinates

### 2. **Homogenization** - Effective Property Calculations

Methods for computing effective properties of heterogeneous materials:

- **Eshelby Tensor** - Analytical solutions for ellipsoidal, cylindrical, and layered inclusions
- **Multi-phase Models** - Cylinder, ellipsoid, and layer assembly homogenization schemes
- **Micromechanical Schemes** - Mori-Tanaka, self-consistent, dilute approximations

### 3. **Material** - Material Characterization

Functions for material property characterization and orientation distribution:

- **ODF** (Orientation Distribution Function) - Crystallographic texture representation
- **PDF** (Probability Density Function) - Statistical phase distributions
- **Crystallography** - Crystal structure and slip system definitions
- **Variant** - Crystallographic variant management for phase transformations

### 4. **Micromechanics** - Multi-scale Material Modeling

Multi-scale modeling frameworks connecting microscopic and macroscopic scales:

- **Phase Characteristics** - Complete material phase definition (geometry, properties, state)
- **Concentration Tensors** - Stress and strain concentration factors
- **Localization** - Microscopic field computation from macroscopic loading

### 5. **UMAT** - Material Constitutive Models

User Material subroutines (UMAT) for finite element analysis, organized by strain formulation and physical behavior:

#### **Mechanical Models** (Small Strain)

**Elasticity:**
- Isotropic, Orthotropic, Transversely Isotropic

**Plasticity:**
- `plastic_isotropic_ccp` - J2 plasticity with isotropic hardening
- `plastic_kin_iso_ccp` - Combined kinematic-isotropic hardening
- `plastic_chaboche_ccp` - Chaboche multi-kinematic hardening model
- `Hill_isoh` - Hill anisotropic plasticity with isotropic hardening
- `Hill_isoh_Nfast` - Hill plasticity with multiple hardening mechanisms
- `Hill_chaboche_ccp` - Hill plasticity with Chaboche hardening
- `Ani_chaboche_ccp` - General anisotropic plasticity with Chaboche hardening
- `DFA_chaboche_ccp` - Distortion-based anisotropy with Chaboche hardening
- `Generic_chaboche_ccp` - Generic anisotropic plasticity framework

**Viscoelasticity:**
- `Zener_fast` - Single Zener (Standard Linear Solid) element
- `Zener_Nfast` - Generalized Maxwell model with N Zener elements
- `Prony_Nfast` - Prony series representation

**Damage:**
- `damage_LLD_0` - Lemaitre ductile damage model
- `damage_weibull` - Weibull statistical damage model

**Shape Memory Alloys (SMA):**
- `SMA_mono` - Monolithic SMA model
- `SMA_mono_cubic` - Cubic crystal SMA model
- `unified_T` - Unified thermomechanical SMA model

**Combined:**
- `Prony_Nfast_Plastic` - Viscoelastoplasticity

**External:**
- `external_umat` - Plugin API for custom material models

#### **Finite Strain Models**

- `neo_hookean_comp` - Compressible Neo-Hookean hyperelasticity
- `neo_hookean_incomp` - Incompressible Neo-Hookean hyperelasticity
- `mooney_rivlin` - Mooney-Rivlin hyperelasticity
- `saint_venant` - Saint-Venant Kirchhoff model
- `generic_hyper_invariants` - Generic hyperelasticity based on invariants
- `generic_hyper_pstretch` - Generic hyperelasticity based on principal stretches
- `hypoelastic_orthotropic` - Hypoelastic orthotropic model

#### **Thermomechanical Models**

Coupled thermomechanical versions of elasticity, plasticity, viscoelasticity, and SMA models.

### 6. **Common Features Across UMATs**

All UMAT implementations share common characteristics:

- **Convex Cutting Plane (CCP) Algorithm** - Robust return mapping for plasticity
- **Consistent Tangent Modulus** - Ensures quadratic convergence in implicit FE
- **Objective Integration** - Frame-invariant stress updates for large rotations
- **Thermal Expansion** - Built-in thermal strain effects
- **Energy Dissipation Tracking** - Mechanical work decomposition

## Mathematical Foundations

### Voigt Notation Convention

Simcoon uses the following Voigt notation for symmetric tensors:

**Stress/Strain Vector:**
\f[
\boldsymbol{\sigma} = [\sigma_{11}, \sigma_{22}, \sigma_{33}, \sigma_{12}, \sigma_{13}, \sigma_{23}]^T
\f]

**Engineering Shear Strains:**

For strain tensors, shear components use engineering notation:
\f[
\boldsymbol{\varepsilon} = [\varepsilon_{11}, \varepsilon_{22}, \varepsilon_{33}, \gamma_{12}, \gamma_{13}, \gamma_{23}]^T
\f]
where \f$ \gamma_{ij} = 2\varepsilon_{ij} \f$ for \f$ i \neq j \f$.

**Stiffness Tensor Convention:**

The elastic stiffness tensor \f$ \mathbf{L} \f$ relates stress and strain as:
\f[
\boldsymbol{\sigma} = \mathbf{L} : \boldsymbol{\varepsilon}
\f]

### Identity Tensors

Several fourth-order identity tensors are provided:

- \f$ \mathbf{I}_{real} \f$ - Standard identity with 0.5 on shear diagonal
- \f$ \mathbf{I}_{vol} \f$ - Volumetric projection tensor
- \f$ \mathbf{I}_{dev} \f$ - Deviatoric projection tensor
- \f$ \widehat{\mathbf{I}} \f$ - Identity with 2 on shear diagonal (for L·M̂=I)

### Stress Measures

Multiple stress measures for finite strain analysis:

- **Cauchy Stress** \f$ \boldsymbol{\sigma} \f$ - True stress in current configuration
- **First Piola-Kirchhoff** \f$ \mathbf{P} \f$ - Two-point tensor
- **Second Piola-Kirchhoff** \f$ \mathbf{S} \f$ - Material stress measure
- **Kirchhoff Stress** \f$ \boldsymbol{\tau} = J\boldsymbol{\sigma} \f$ - Weighted Cauchy stress
- **Biot Stress** - Work-conjugate to Biot strain

### Strain Measures

Various strain measures for different kinematic descriptions:

- **Green-Lagrange** \f$ \mathbf{E} = \frac{1}{2}(\mathbf{C} - \mathbf{I}) \f$
- **Euler-Almansi** \f$ \mathbf{e} = \frac{1}{2}(\mathbf{I} - \mathbf{b}^{-1}) \f$
- **Logarithmic (Hencky)** \f$ \mathbf{H} = \ln(\mathbf{V}) \f$
- **Engineering Strain** - Small strain approximation

## Usage Patterns

### Basic Elastic Stiffness Construction

```cpp
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>

using namespace simcoon;

// Isotropic material
double E = 70000;   // MPa
double nu = 0.3;
mat L = L_iso(E, nu, "Enu");

// Orthotropic material
double E1 = 150000, E2 = 10000, E3 = 10000;
double nu12 = 0.25, nu13 = 0.25, nu23 = 0.4;
double G12 = 5000, G13 = 5000, G23 = 3500;
mat L_ortho = L_ortho(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, "EnuG");
```

### Stress Invariants and Flow Directions

```cpp
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>

using namespace simcoon;

vec sigma = {100, 50, 30, 10, 5, 0};  // Stress state

// von Mises equivalent stress
double sigma_eq = Mises_stress(sigma);

// Flow direction for J2 plasticity
vec eta = eta_stress(sigma);  // Returns 3/2 * s / sigma_eq
```

### UMAT Integration

```cpp
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>

using namespace simcoon;

// Material properties: [E, nu, alpha, sigma_Y, k, m]
vec props = {70000, 0.3, 1e-5, 200, 500, 0.2};

// State variables: [T_init, p, EP_11, EP_22, EP_33, EP_12, EP_13, EP_23]
vec statev = zeros(8);

// Strain increment
vec DEtot = {0.001, 0, 0, 0, 0, 0};

// Call UMAT
umat_plasticity_iso_CCP(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                        6, props, 8, statev, T, DT, Time, DTime,
                        Wm, Wm_r, Wm_ir, Wm_d, 3, 3, false, 0, tnew_dt);
```

## Design Philosophy

### Modularity

Functions are organized by mathematical operation rather than specific material models, promoting reusability across different constitutive laws.

### Performance

- Armadillo library for optimized linear algebra
- Voigt notation for efficient 6×6 matrix operations
- Minimal dynamic allocation in performance-critical paths

### Robustness

- CCP algorithm for plasticity provides superior convergence
- Fischer-Burmeister complementarity ensures non-negative plastic multipliers
- Consistent tangent moduli for quadratic FE convergence

### Extensibility

- Plugin API for external material models
- Template-based design allows easy addition of new models
- Separation of kinematics, constitutive laws, and solution algorithms

## References

1. **Convex Cutting Plane Algorithm:**
   - Ortiz, M., & Simo, J. C. (1986). "An analysis of a new class of integration algorithms for elastoplastic constitutive relations." *International Journal for Numerical Methods in Engineering*.

2. **Consistent Tangent Modulus:**
   - Simo, J. C., & Taylor, R. L. (1985). "Consistent tangent operators for rate-independent elastoplasticity." *Computer Methods in Applied Mechanics and Engineering*.

3. **Hyperelasticity:**
   - Holzapfel, G. A. (2000). *Nonlinear Solid Mechanics: A Continuum Approach for Engineering*. Wiley.

4. **Micromechanics:**
   - Mori, T., & Tanaka, K. (1973). "Average stress in matrix and average elastic energy of materials with misfitting inclusions." *Acta Metallurgica*.
   - Eshelby, J. D. (1957). "The determination of the elastic field of an ellipsoidal inclusion." *Proceedings of the Royal Society A*.

## See Also

- [Simulation Module](simulation_overview.md) - Solver framework and parameter identification
- [Function Reference](functions/) - Detailed API documentation
- [UMAT Reference](umat/) - Material model implementations
