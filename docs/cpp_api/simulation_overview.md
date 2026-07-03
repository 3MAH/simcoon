# Simulation Module Overview

## Introduction

The Simulation module provides the infrastructure for executing multi-physics simulations. It includes the solver framework, phase management, and mathematical utilities.

## Module Organization

### 1. **Solver** - Simulation Execution Framework

The solver manages the execution of finite element-like simulations with support for mechanical and thermomechanical loading paths.

#### **Key Components:**

- **solver.hpp** - Main simulation driver
- **step.hpp** - Base class for loading steps
- **step_meca.hpp** - Mechanical loading steps
- **step_thermomeca.hpp** - Thermomechanical coupled loading steps
- **block.hpp** - Block structure for multi-steop simulations, incuding repetitive cycles
- **output.hpp** - Results output management
- **read.hpp** - Input file parsing

#### **Solver Capabilities:**

- **Loading Conditions:**
  - Strain control
  - Stress control
  - Mixed stress-strain control
  - Temperature-prescribed paths
  - Coupled thermomechanical loading

- **Solution Methods:**
  - Implicit Newton-Raphson
  - Adaptive time stepping

- **Output Control:**
  - Stress-strain results using various measures
  - Internal state variables evolution
  - Energy dissipation tracking
  - Custom output fields

### 2. **Phase** - Material Phase Management

Multi-scale material representation with support for heterogeneous microstructures.

#### **Core Classes:**

##### **phase_characteristics**
Complete representation of a material phase including:
- Geometry (shape, size, orientation)
- Material properties (UMAT name, parameters)
- State variables (stress, strain, internal variables)
- Concentration tensors (localization operators)
- Subphases for hierarchical multi-scale modeling

##### **state_variables**
Container for mechanical state at a material point:
- Strain measures (engineering, Green-Lagrange, logarithmic)
- Stress measures (Cauchy, Kirchhoff, Piola-Kirchhoff)
- Deformation gradient, stretch tensors, rotation
- Internal state variables (UMAT-specific)
- Temperature field

**Variants:**
- **state_variables_M** - Pure mechanical state
- **state_variables_T** - Thermomechanical coupled state

##### **material_characteristics**
Material property container:
- UMAT identifier
- Material parameters (props vector)
- Number of state variables required
- Material orientation (Euler angles)

##### **geometry**
Base class for inclusion shapes:
- **ellipsoid** - Ellipsoidal inclusions (fibers, particles)
- **cylinder** - Cylindrical inclusions (infinite fibers)
- **layer** - Layered structures (laminates)

##### **output**
Results management:
- Selected output fields
- Output frequency control
- File format selection

##### **I/O Operations:**
- **read.hpp** - Phase file parsing
- **write.hpp** - State serialization

### 3. **Maths** - Mathematical Utilities

Mathematical tools supporting simulation:

#### **rotation.hpp**
Rotation operations for objective stress integration:
- Rotation matrices from Euler angles
- Rotation of vectors and tensors
- Quaternion operations
- Active/passive rotation conventions

**Key Functions:**
- `Rotation::from_euler()` - Construct rotation from Euler angles (ZXZ, ZYZ, XYZ, etc.)
- `Rotation::from_axis_angle()` - Rotation from axis-angle representation
- `rotate_strain()` - Rotate strain tensor (Voigt notation)
- `rotate_stress()` - Rotate stress tensor (Voigt notation)
- `rotate_stiffness()` - Rotate stiffness matrix
- `rotate_compliance()` - Rotate compliance matrix

#### **lagrange.hpp**
Lagrange multiplier methods for constrained problems:
- Exponential penalty functions
- Power-law penalties
- Derivatives for numerical implementations

**Applications:**
- Physical constraints (volumic fraction of phases)

#### **stats.hpp**
Statistical functions for data analysis:
- Mean, variance, standard deviation
- Probability distributions

#### **solver.hpp**
Linear system solvers:
- Quadratic functions

#### **num_solve.hpp**
Nonlinear equation solvers:
- Newton-Raphson method
- Modified Newton methods using Fischer-Burmeister complementarity functions

**Fischer-Burmeister Function:**

For complementarity problems \f$ a \geq 0, b \geq 0, a \cdot b = 0 \f$:
\f[
FB(a,b) = \sqrt{a^2 + b^2} - a - b = 0
\f]

Used extensively in plasticity return mapping.

#### **random.hpp**
Random number generation:
- Uniform distributions
- Normal (Gaussian) distributions
- Seeding control
- Reproducibility support

### 4. **Geometry** - Geometric Primitives

Geometric representations for multi-phase materials:

#### **geometry (base class)**
Abstract interface for inclusion shapes:
```cpp
class geometry {
  public:
    int shape_type;     // 1=ellipsoid, 2=cylinder, 3=layer
    vec params;         // Shape parameters
    double concentration;  // Volume fraction
};
```

#### **ellipsoid**
Ellipsoidal inclusions characterized by semi-axes \f$ (a_1, a_2, a_3) \f$:
- Aspect ratios define shape (sphere, oblate, prolate)
- Orientation defined by Euler angles
- Eshelby tensor computation

#### **cylinder**
Infinitely long cylindrical fibers:
- Characterized by aspect ratio \f$ a_1/a_2 \f$ of cross-section
- Fiber direction alignment
- Transverse isotropy

#### **layer**
Planar layered structures:
- Layer thickness
- Normal direction
- Perfect bonding assumption

## Simulation Workflow

### 1. Define Material Properties

```cpp
// Define a phase
phase_characteristics phase1;
phase1.sptr_matprops->umat_name = "ELISO";
phase1.sptr_matprops->props = {70000, 0.3, 1e-5};  // E, nu, alpha
phase1.sptr_shape->shape_type = 1;  // Ellipsoid
phase1.concentration = 0.6;
```

### 2. Create Loading Path

```cpp
// Create mechanical loading step
step_meca step1;
step1.number = 1;
step1.Dn_init = 100;  // Number of increments
step1.Dn_mini = 10;
step1.Dn_maxi = 1000;

// Define strain control in direction 11
step1.cBC_meca = "E";
step1.BC_meca(0) = 0.01;  // 1% strain
```

### 3. Execute Simulation

```cpp
solver(umat_name, props, nstatev, psi, theta, phi);
```

### 4. Post-process Results

Results are written to output files specified in the control file.

## Advanced Features

### Multi-scale Modeling

Hierarchical material definition with nested phases:

```cpp
phase_characteristics matrix, inclusion, composite;

// Define inclusion within matrix
matrix.sub_phases.push_back(inclusion);

// Define composite with matrix containing inclusions
composite.sub_phases.push_back(matrix);
```

### Adaptive Time Stepping

The solver automatically adjusts time step size based on:
- Convergence rate
- Error estimates
- User-defined bounds

### Parallelization

Some operations support parallel execution (using OpenMP):
- Multiple phase response

## Configuration Files

### Material Properties File

Format:
```
Number_of_phases
Phase_1:
  UMAT_name
  Number_of_properties
  prop_1 prop_2 ... prop_n
  Number_of_statev
Phase_2:
  ...
```

### Path File

Defines loading sequence:
```
Number_of_steps
Step_1:
  Type (mechanical/thermomechanical)
  Number_of_increments
  Control_type (E/S/T)
  BC_1 BC_2 BC_3 BC_4 BC_5 BC_6
  Temperature Temperature_increment
Step_2:
  ...
```

### Output File

Configures results output:
```
Number_of_output_blocks
Block_1:
  Output_type (stress/strain/STATEV)
  Components (space-separated indices)
  Frequency
```

## Performance Considerations

### Memory Management

- State variables allocated once per increment
- Minimal dynamic allocation in tight loops
- Smart pointers for phase hierarchies

### Computational Efficiency

- Efficient Voigt notation operations (6×6 vs 3×3×3×3)
- Armadillo's optimized BLAS/LAPACK backend
- Sparse matrix support for large systems

### Numerical Stability

- Consistent tangent moduli prevent drift
- Adaptive time stepping for stiff problems
- Fischer-Burmeister for robust complementarity

## References

1. **Micromechanics and Multi-scale Modeling:**
   - Qu, J., & Cherkaoui, M. (2006). *Fundamentals of Micromechanics of Solids*. Wiley.

2. **Fischer-Burmeister Method:**
   - Fischer, A. (1997). "Solution of monotone complementarity problems with locally Lipschitzian functions." *Mathematical Programming*.

## See Also

- [Continuum Mechanics Module](continuum_mechanics_overview.md) - Material models and constitutive functions
- [Solver API](solver/) - Detailed solver documentation
