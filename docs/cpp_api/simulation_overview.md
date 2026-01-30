# Simulation Module Overview

## Introduction

The Simulation module provides the infrastructure for executing multi-physics simulations, parameter identification, and optimization. It includes the solver framework, phase management, mathematical utilities, and identification algorithms for material parameter calibration.

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
- **read_json.hpp** - JSON-based phase file parsing
- **write.hpp** - State serialization

Phase configurations use JSON format. See the
[Micromechanics I/O documentation](../simulation/micromechanics.rst) for the Python API
and JSON format specifications.

### 3. **Identification** - Parameter Identification Framework

Inverse analysis tools for material parameter calibration from experimental data.

#### **Optimization Framework:**

##### **identification.hpp**
Main identification driver supporting multiple optimization algorithms:
- Genetic Algorithms (GA)
- Gradient-based methods (Levenberg-Marquardt)
- Hybrid strategies

##### **Core Components:**

**parameters.hpp**
Defines optimization variables:
```cpp
class parameters {
  public:
    double value;           // Current value
    double min_value;       // Lower bound
    double max_value;       // Upper bound
    string key;             // Parameter identifier
    int ninit;              // Number of initializations
};
```

**constants.hpp**
Fixed values during optimization:
```cpp
class constants {
  public:
    double value;
    string key;
    int ninit;
};
```

**individual.hpp**
Solution candidate in population-based methods:
- Parameter vector
- Cost function value
- Constraint violations
- Fitness ranking

**generation.hpp**
Population management for evolutionary algorithms:
- Population initialization
- Selection operators
- Crossover and mutation
- Elitism strategies

**methods.hpp**
Optimization algorithm implementations:
- Cost function evaluation
- Gradient computation (numerical/analytical)
- Hessian approximation

**optimize.hpp**
High-level optimization loop:
- Convergence checking
- Iteration management
- History tracking
- Checkpointing

**opti_data.hpp**
Experimental data management:
- Loading experimental files
- Data interpolation
- Weight assignment

**doe.hpp** (Design of Experiments)
Sampling strategies for parameter space exploration:
- Latin Hypercube Sampling (LHS)
- Random sampling

**script.hpp**
Script interpretation for identification workflows.

#### **Cost Function Definition:**

The identification minimizes a weighted sum of squared residuals:

\f[
f(\mathbf{p}) = \sum_{i=1}^{N_{exp}} \sum_{j=1}^{N_{pts}} w_{ij} \left( y^{exp}_{ij} - y^{sim}_{ij}(\mathbf{p}) \right)^2
\f]

where:
- \f$ \mathbf{p} \f$ = parameter vector
- \f$ N_{exp} \f$ = number of experimental datasets
- \f$ N_{pts} \f$ = number of data points per dataset
- \f$ w_{ij} \f$ = weight for point j in experiment i
- \f$ y^{exp}_{ij} \f$ = experimental observation
- \f$ y^{sim}_{ij}(\mathbf{p}) \f$ = simulation prediction

### 4. **Maths** - Mathematical Utilities

Mathematical tools supporting simulation and identification:

#### **rotation.hpp**
Rotation operations for objective stress integration:
- Rotation matrices from Euler angles
- Rotation of vectors and tensors
- Quaternion operations
- Active/passive rotation conventions

**Key Functions:**
- `fillR_euler()` - Construct rotation matrix from Euler angles (ZXZ, ZYZ conventions)
- `fillR_angle()` - Rotation matrix from axis-angle representation
- `rotate_strain()` - Rotate strain tensor (Voigt notation)
- `rotate_stress()` - Rotate stress tensor (Voigt notation)
- `rotateL()` - Rotate stiffness matrix
- `rotateM()` - Rotate compliance matrix

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

### 5. **Geometry** - Geometric Primitives

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

**Note:** The C++ `solver()` function has been replaced by the Python Solver API.
Use `simcoon.solver.Solver` in Python for all new simulations. See the
[Python Solver documentation](../simulation/solver.rst) for details.

### 4. Post-process Results

Results are returned as a list of `StateVariablesM` objects in Python.

## Parameter Identification Workflow

**Note:** The C++ identification module is deprecated in v2.0.
Use Python with `scipy.optimize` for parameter identification:

```python
from scipy.optimize import minimize, least_squares
from simcoon.solver import Solver, Block, StepMeca
import numpy as np

def simulate(params):
    """Run simulation with given parameters."""
    E, sigma_Y = params
    props = np.array([E, 0.3, sigma_Y, ...])

    step = StepMeca(DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]))
    block = Block(steps=[step], umat_name='EPICP', props=props, nstatev=8)
    solver = Solver(blocks=[block])

    history = solver.solve()
    return np.array([h.sigma[0] for h in history])

def cost_function(params, exp_data):
    """Compute cost between simulation and experiment."""
    sim_data = simulate(params)
    return np.sum((sim_data - exp_data)**2)

# Run optimization
x0 = np.array([200000, 300])  # Initial guess: E, sigma_Y
bounds = [(50000, 300000), (100, 500)]
result = minimize(cost_function, x0, args=(exp_data,), bounds=bounds)
print(f"Optimal parameters: E={result.x[0]:.0f}, sigma_Y={result.x[1]:.1f}")
```

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
- Population-based optimization (outdated)

### Checkpointing

Identification can be checkpointed for:
- Recovery from interruptions
- Continuation of long-running optimizations
- Sensitivity analysis

## Configuration Files

**Note:** The file-based configuration format is deprecated in v2.0.
Use JSON configuration with the Python API instead. See [Solver Documentation](../simulation/solver.rst).

### JSON Configuration (Recommended)

Material configuration (`material.json`):
```json
{
  "name": "ELISO",
  "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5},
  "nstatev": 1,
  "orientation": {"psi": 0, "theta": 0, "phi": 0}
}
```

Path configuration (`path.json`):
```json
{
  "initial_temperature": 293.15,
  "blocks": [
    {
      "type": "mechanical",
      "control_type": "small_strain",
      "ncycle": 1,
      "steps": [
        {
          "time": 30.0,
          "Dn_init": 1.0,
          "Dn_inc": 0.01,
          "DEtot": [0.01, 0, 0, 0, 0, 0],
          "Dsigma": [0, 0, 0, 0, 0, 0],
          "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
          "DT": 0
        }
      ]
    }
  ]
}
```

### Legacy File Formats (Deprecated)

The following legacy file formats are deprecated:

- `material.dat` - Text-based material properties
- `path.txt` - Text-based loading path
- `output.dat` - Text-based output configuration
- `solver_essentials.inp` - Solver type configuration
- `solver_control.inp` - Solver convergence parameters

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

1. **Genetic Algorithms:**
   - Goldberg, D. E. (1989). *Genetic Algorithms in Search, Optimization, and Machine Learning*. Addison-Wesley.

2. **Micromechanics and Multi-scale Modeling:**
   - Qu, J., & Cherkaoui, M. (2006). *Fundamentals of Micromechanics of Solids*. Wiley.

3. **Parameter Identification:**
   - Mahnken, R., & Stein, E. (1996). "Parameter identification for viscoplastic models based on analytical derivatives of a least-squares functional and stability investigations." *International Journal of Plasticity*.

4. **Fischer-Burmeister Method:**
   - Fischer, A. (1997). "Solution of monotone complementarity problems with locally Lipschitzian functions." *Mathematical Programming*.

## See Also

- [Continuum Mechanics Module](continuum_mechanics_overview.md) - Material models and constitutive functions
- [Solver API](solver/) - Detailed solver documentation
- [Identification API](identification/) - Parameter identification reference
