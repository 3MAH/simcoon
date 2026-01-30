# Simcoon v2.0 Migration Guide

This guide helps you transition from the legacy file-based C++ solver to the new Python-based solver API introduced in Simcoon v2.0.

## Overview of Changes

### What's New in v2.0

- **Python Solver API**: Full-featured Python solver with `Solver`, `Block`, `StepMeca`, and `StepThermomeca` classes
- **JSON Configuration**: Load/save material properties and loading paths as JSON
- **Python Identification Module**: Material parameter identification using scipy.optimize
- **Improved History Access**: Direct access to all state variables at each increment

### What's Deprecated/Removed

| Legacy Function | Status | Replacement |
|-----------------|--------|-------------|
| `sim.solver()` | Removed | `simcoon.solver.Solver` class |
| `sim.read_matprops()` | Removed | `simcoon.solver.load_material_json()` |
| `sim.read_path()` | Removed | `simcoon.solver.load_path_json()` |
| `sim.identification()` | Removed | `simcoon.identification` module |
| `sim.calc_cost()` | Removed | `simcoon.identification.mse()`, etc. |

## Migration Examples

### Basic Simulation

**Legacy (v1.x):**
```python
import simcoon as sim
import numpy as np

umat_name = "ELISO"
props = np.array([210000.0, 0.3, 1e-5])
nstatev = 1

# Requires external path.txt file and output directory
sim.solver(
    umat_name, props, nstatev,
    0.0, 0.0, 0.0,  # Euler angles
    0,  # solver_type
    2,  # corate_type
    "data", "results",
    "path.txt", "output.txt"
)

# Read results from file
data = np.loadtxt("results/output_global-0.txt")
```

**New (v2.0):**
```python
import numpy as np
from simcoon.solver import Solver, Block, StepMeca

props = np.array([210000.0, 0.3, 1e-5])

# Define loading directly in Python
step = StepMeca(
    DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=50
)

block = Block(
    steps=[step],
    umat_name="ELISO",
    props=props,
    nstatev=1,
    control_type='small_strain',
    corate_type='logarithmic'
)

solver = Solver(blocks=[block])
history = solver.solve()

# Access results directly
strain = np.array([h.Etot[0] for h in history])
stress = np.array([h.sigma[0] for h in history])
```

### Multi-Step Loading

**Legacy (v1.x):**
Required editing `path.txt` file manually.

**New (v2.0):**
```python
# Tension
step1 = StepMeca(
    DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=100
)

# Compression
step2 = StepMeca(
    DEtot_end=np.array([-0.04, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=200
)

block = Block(
    steps=[step1, step2],
    umat_name="EPICP",
    props=props,
    nstatev=8
)

solver = Solver(blocks=[block])
history = solver.solve()
```

### JSON Configuration

**New (v2.0):**
```python
from simcoon.solver import load_material_json, load_path_json, Solver

# Load from JSON files
material = load_material_json("material.json")
blocks = load_path_json("path.json")

# Update blocks with material
for block in blocks:
    block.umat_name = material['umat_name']
    block.props = np.array(material['props'])
    block.nstatev = material['nstatev']

solver = Solver(blocks=blocks)
history = solver.solve()
```

### Material Parameter Identification

**Legacy (v1.x):**
```python
# Required specific file structure and C++ bindings
sim.identification()  # Limited Python control
```

**New (v2.0):**
```python
from simcoon.identification import (
    IdentificationProblem,
    levenberg_marquardt,
    differential_evolution
)
from simcoon.solver import Solver, Block, StepMeca

# Define simulation function
def simulate(params):
    E, sigma_Y = params
    props = np.array([E, 0.3, 0.0, sigma_Y, 1000, 0.3])

    step = StepMeca(DEtot_end=np.array([0.05, 0, 0, 0, 0, 0]))
    block = Block(steps=[step], umat_name="EPICP", props=props, nstatev=8)
    solver = Solver(blocks=[block])
    history = solver.solve()

    return {'stress': np.array([h.sigma[0] for h in history])}

# Create identification problem
problem = IdentificationProblem(
    parameters=[
        {'name': 'E', 'bounds': (150000, 250000)},
        {'name': 'sigma_Y', 'bounds': (200, 500)},
    ],
    simulate=simulate,
    exp_data={'stress': experimental_stress_data},
)

# Run optimization
result = levenberg_marquardt(problem)
print(f"Identified E: {result.x[0]:.0f}")
print(f"Identified sigma_Y: {result.x[1]:.0f}")
```

## Control Type and Corate Type Mappings

**Control Types:**
| String | Integer | Description |
|--------|---------|-------------|
| `'small_strain'` | 1 | Small strain formulation |
| `'green_lagrange'` | 2 | Green-Lagrange strain |
| `'logarithmic'` | 3 | Logarithmic (Hencky) strain |
| `'biot'` | 4 | Biot strain |
| `'F'` | 5 | Deformation gradient control |
| `'gradU'` | 6 | Displacement gradient control |

**Corate Types:**
| String | Integer | Description |
|--------|---------|-------------|
| `'jaumann'` | 0 | Jaumann (Zaremba-Jaumann) rate |
| `'green_naghdi'` | 1 | Green-Naghdi rate |
| `'logarithmic'` | 2 | Logarithmic rate |
| `'logarithmic_R'` | 3 | Logarithmic rate with rotation |
| `'truesdell'` | 4 | Truesdell rate |
| `'logarithmic_F'` | 5 | Logarithmic rate from F |

## State Variables Access

**Legacy (v1.x):**
Results written to files with fixed column format.

**New (v2.0):**
```python
history = solver.solve()

# Each entry in history is a StateVariablesM object
for state in history:
    # Strains (Voigt notation)
    print(state.Etot)    # Total strain
    print(state.DEtot)   # Strain increment

    # Stresses
    print(state.sigma)   # Cauchy stress
    print(state.PKII)    # 2nd Piola-Kirchhoff stress

    # Deformation
    print(state.F0)      # Deformation gradient (start)
    print(state.F1)      # Deformation gradient (end)
    print(state.R)       # Rotation tensor

    # Work quantities
    print(state.Wm)      # [Wm, Wm_r, Wm_ir, Wm_d]

    # Internal variables
    print(state.statev)

    # Tangent stiffness
    print(state.Lt)      # 6x6 tangent modulus
```

## Micromechanics (Unchanged)

The homogenization functions remain unchanged:

```python
import simcoon as sim

# These still work the same way
L_eff = sim.L_eff(umat_name, props, nstatev, psi, theta, phi)
S = sim.Eshelby_sphere(nu)
S = sim.Eshelby_prolate(nu, aspect_ratio)
```

## Common Issues

### "AttributeError: module 'simcoon' has no attribute 'solver'"

This error occurs when calling the old `sim.solver()` function. Use the new API:
```python
from simcoon.solver import Solver, Block, StepMeca
```

### "TypeError: solver() missing required argument"

The new `Solver` class doesn't take the same arguments. See examples above.

### Import Errors for Identification

If you get import errors for `identification`, make sure you have scipy installed:
```bash
pip install scipy
```

## Need Help?

- [Documentation](https://3mah.github.io/simcoon-docs/)
- [Examples](https://github.com/3MAH/simcoon/tree/master/examples)
- [Issues](https://github.com/3MAH/simcoon/issues)
