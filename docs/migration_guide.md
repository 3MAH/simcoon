# Simcoon v2.0 Migration Guide

A comprehensive guide for transitioning from Simcoon v1.x to v2.0.

---

## Executive Summary

Simcoon v2.0 introduces significant API changes focused on improving usability and flexibility:

| Area | v1.x Approach | v2.0 Approach |
|------|---------------|---------------|
| **Solver** | `sim.solver()` with file-based I/O | `simcoon.solver.Solver` class with Python objects |
| **Configuration** | Text files (`path.txt`, `material.dat`) | JSON files or Python dataclasses |
| **Identification** | C++ `identification()` with file workflow | Pure Python `simcoon.identification` module |
| **Micromechanics** | `Nphases.dat`, `Nlayers.dat` text files | Python dataclasses with JSON I/O |
| **Results** | Output files with fixed column format | Direct access to `HistoryPoint` objects |

### Key Benefits of v2.0

- **Full Python control**: Define simulations programmatically without external files
- **Better debugging**: Step through solver iterations in Python debugger
- **Flexible I/O**: JSON format for human-readable, version-controllable configurations
- **Modern optimization**: scipy/sklearn integration for parameter identification
- **Type safety**: Python dataclasses with clear type hints

---

## Table of Contents

1. [Solver API Migration](#1-solver-api-migration)
2. [File-based to JSON-based Configuration](#2-file-based-to-json-based-configuration)
3. [Python Identification Module](#3-python-identification-module)
4. [Micromechanics Data Classes](#4-micromechanics-data-classes)
5. [API Comparison Tables](#5-api-comparison-tables)
6. [Common Migration Patterns](#6-common-migration-patterns)
7. [Troubleshooting](#7-troubleshooting)

---

## 1. Solver API Migration

### v1.x: File-based Solver

The legacy solver required external configuration files and wrote results to output files:

```python
# v1.x (DEPRECATED)
import simcoon as sim
import numpy as np

umat_name = "ELISO"
props = np.array([210000.0, 0.3, 1e-5])
nstatev = 1

# Required: data/path.txt file
# Required: output directory
sim.solver(
    umat_name, props, nstatev,
    0.0, 0.0, 0.0,      # Euler angles (psi, theta, phi)
    0,                   # solver_type
    2,                   # corate_type (integer)
    "data", "results",   # input/output directories
    "path.txt", "output.txt"
)

# Results written to files - must parse manually
data = np.loadtxt("results/output_global-0.txt")
strain = data[:, 0]
stress = data[:, 1]
```

### v2.0: Python Object-Oriented Solver

The new API uses Python classes for complete programmatic control:

```python
# v2.0 (RECOMMENDED)
import numpy as np
from simcoon.solver import Solver, Block, StepMeca

props = np.array([210000.0, 0.3, 1e-5])

# Define loading step
step = StepMeca(
    DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=50,
    Dn_mini=1,
    Dn_inc=100,
    time=1.0
)

# Create block with material
block = Block(
    steps=[step],
    umat_name="ELISO",
    props=props,
    nstatev=1,
    control_type='small_strain',  # String instead of integer
    corate_type='logarithmic'     # String instead of integer
)

# Run solver
solver = Solver(blocks=[block])
history = solver.solve()

# Direct access to results
strain = np.array([h.Etot[0] for h in history])
stress = np.array([h.sigma[0] for h in history])
```

### Key Class Reference

#### `Solver` Class

```python
from simcoon.solver import Solver

solver = Solver(
    blocks=[block1, block2],  # List of Block objects
    max_iter=10,              # Newton-Raphson max iterations
    tol=1e-9,                 # Convergence tolerance
    lambda_solver=10000.0     # Stiffness for strain-controlled components
)

history = solver.solve(sv_init=None)  # Optional initial state
```

#### `Block` Class

```python
from simcoon.solver import Block

block = Block(
    steps=[step1, step2],     # List of StepMeca/StepThermomeca
    nstatev=1,                # Number of state variables
    umat_name="ELISO",        # UMAT name (5 chars)
    umat_type="mechanical",   # "mechanical" or "thermomechanical"
    props=np.array([...]),    # Material properties
    control_type='small_strain',  # See control types below
    corate_type='jaumann',        # See corate types below
    ncycle=1                  # Number of cycles to repeat
)
```

#### `StepMeca` Class

```python
from simcoon.solver import StepMeca

step = StepMeca(
    DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),  # Target strain increment
    Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),    # Target stress increment
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    time=1.0,        # Step duration
    Dn_init=1,       # Initial increment count
    Dn_mini=1,       # Minimum increments
    Dn_inc=100       # Maximum increments
)
```

#### `StepThermomeca` Class

```python
from simcoon.solver import StepThermomeca

step = StepThermomeca(
    DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    DT_end=50.0,                    # Temperature increment
    Q_end=0.0,                      # Heat flux
    thermal_control='temperature'   # 'temperature' or 'heat_flux'
)
```

---

## 2. File-based to JSON-based Configuration

### v1.x: Text File Format

**path.txt** (tab-separated, fixed columns):
```
#Initial_temperature
293.15
#Number_of_blocks
1
#Block
1	1	0	0	2
#Number_of_steps
1
#Steps
2	0	30	1	1	0.01	0.01	0	0	0	0	0	0	1	0	0	0	0	0	0
```

**material.dat** (tab-separated):
```
ELISO	1	210000	0.3	1e-5	0	0	0
```

**Nphases.dat** (tab-separated):
```
Number	Name	save	c	psi	theta	phi	nstatev	nprops	props
0	ELISO	1	0.7	0	0	0	1	3	3500	0.35	6e-5
1	ELISO	1	0.3	0	0	0	1	3	72000	0.22	5e-6
```

### v2.0: JSON Format

**material.json**:
```json
{
  "name": "ELISO",
  "props": {"E": 210000, "nu": 0.3, "alpha": 1e-5},
  "nstatev": 1,
  "orientation": {"psi": 0, "theta": 0, "phi": 0}
}
```

**path.json**:
```json
{
  "initial_temperature": 293.15,
  "blocks": [
    {
      "type": "mechanical",
      "control_type": "small_strain",
      "corate_type": "logarithmic",
      "ncycle": 1,
      "steps": [
        {
          "time": 30.0,
          "Dn_init": 1,
          "Dn_mini": 1,
          "Dn_inc": 100,
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

### JSON I/O Functions

```python
from simcoon.solver.io import (
    load_material_json, save_material_json,
    load_path_json, save_path_json,
    load_simulation_json  # Combined loading
)

# Load material
material = load_material_json('material.json')
# Returns: {'name': 'ELISO', 'props': array([...]), 'nstatev': 1, 'orientation': {...}}

# Load path (returns Block objects)
path = load_path_json('path.json')
# Returns: {'initial_temperature': 293.15, 'blocks': [Block(...), ...]}

# Combined loading (assigns material to blocks)
sim = load_simulation_json('material.json', 'path.json')
solver = Solver(blocks=sim['blocks'])
history = solver.solve()

# Save configurations
save_material_json('material.json', 'ELISO', props, nstatev=1)
save_path_json('path.json', blocks, initial_temperature=293.15)
```

---

## 3. Python Identification Module

### v1.x: C++ Identification (DEPRECATED)

```python
# v1.x (DEPRECATED)
import simcoon as sim

# Required specific file structure:
# - data/id_params.txt
# - data/exp_data.txt
# - data/path.txt

sim.identification()  # Limited Python control, file-based workflow
```

### v2.0: Pure Python Module

```python
# v2.0 (RECOMMENDED)
from simcoon.identification import (
    IdentificationProblem,
    ParameterSpec,
    OptimizationResult,
    levenberg_marquardt,
    differential_evolution,
    hybrid_optimization,
    nelder_mead,
    mse, mae, r2, weighted_mse, huber_loss,
    compute_sensitivity,
    compute_jacobian,
    correlation_matrix,
)
from simcoon.solver import Solver, Block, StepMeca
import numpy as np
```

### Complete Identification Example

```python
import numpy as np
from simcoon.identification import IdentificationProblem, levenberg_marquardt
from simcoon.solver import Solver, Block, StepMeca

# Experimental data (strain-stress curve)
exp_strain = np.linspace(0, 0.05, 100)
exp_stress = np.array([...])  # Your experimental data

# Define simulation function
def simulate(params):
    E, sigma_Y, H = params
    props = np.array([E, 0.3, 0.0, sigma_Y, H, 0.3])

    step = StepMeca(
        DEtot_end=np.array([0.05, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=100
    )
    block = Block(steps=[step], umat_name="EPICP", props=props, nstatev=8)

    solver = Solver(blocks=[block])
    history = solver.solve()

    return {
        'stress': np.array([h.sigma[0] for h in history]),
        'strain': np.array([h.Etot[0] for h in history])
    }

# Define parameters with bounds
problem = IdentificationProblem(
    parameters=[
        {'name': 'E', 'bounds': (150000, 250000), 'initial': 200000},
        {'name': 'sigma_Y', 'bounds': (200, 500), 'initial': 350},
        {'name': 'H', 'bounds': (500, 2000), 'initial': 1000},
    ],
    simulate=simulate,
    exp_data={'stress': exp_stress},
    weights={'stress': 1.0},
    cost_type='mse'
)

# Run Levenberg-Marquardt optimization
result = levenberg_marquardt(problem, verbose=2)

print(f"Identified parameters:")
for name, val in zip(result.parameter_names, result.x):
    print(f"  {name}: {val:.2f}")
print(f"Final cost: {result.cost:.6e}")
print(f"Converged: {result.success}")
```

### Available Optimizers

| Optimizer | Best For | Example |
|-----------|----------|---------|
| `levenberg_marquardt` | Well-posed problems with good initial guess | `levenberg_marquardt(problem, ftol=1e-8)` |
| `differential_evolution` | Global search, many local minima | `differential_evolution(problem, maxiter=500)` |
| `nelder_mead` | Smooth problems, few parameters | `nelder_mead(problem, adaptive=True)` |
| `hybrid_optimization` | Robust global + local search | `hybrid_optimization(problem, n_restarts=5)` |

### Parameter Specification

```python
from simcoon.identification import ParameterSpec

# Using ParameterSpec class
param = ParameterSpec(
    name='E',
    bounds=(100000, 300000),
    initial=200000,        # Optional: defaults to midpoint
    scale=None,            # Optional: defaults to range
    fixed=False            # Set True to hold fixed
)

# Or use dict shorthand
param_dict = {'name': 'E', 'bounds': (100000, 300000), 'initial': 200000}
```

---

## 4. Micromechanics Data Classes

### v1.x: Text File Format

**Nellipsoids.dat**:
```
Number	coatingof	Name	save	c	psi_mat	theta_mat	phi_mat	a1	a2	a3	psi_geo	theta_geo	phi_geo	nstatev	nprops	props
0	0	ELISO	1	0.7	0	0	0	1	1	1	0	0	0	1	3	3500	0.35	6e-5
1	0	ELISO	1	0.3	0	0	0	20	1	1	0	0	0	1	3	72000	0.22	5e-6
```

### v2.0: Python Dataclasses

```python
from simcoon.solver.micromechanics import (
    Phase, Layer, Ellipsoid, Cylinder, Section,
    MaterialOrientation, GeometryOrientation,
    load_phases_json, save_phases_json,
    load_layers_json, save_layers_json,
    load_ellipsoids_json, save_ellipsoids_json,
    load_cylinders_json, save_cylinders_json,
)
import numpy as np
```

### Ellipsoid Example

```python
from simcoon.solver.micromechanics import Ellipsoid, save_ellipsoids_json

# Create matrix phase (spherical)
matrix = Ellipsoid(
    number=0,
    concentration=0.7,
    umat_name="ELISO",
    props=np.array([3500, 0.35, 6e-5]),
    a1=1, a2=1, a3=1,  # Spherical
    nstatev=1
)

# Create fiber phase (prolate spheroid)
fiber = Ellipsoid(
    number=1,
    concentration=0.3,
    umat_name="ELISO",
    props=np.array([72000, 0.22, 5e-6]),
    a1=20, a2=1, a3=1,  # Aspect ratio 20
    geometry_orientation=GeometryOrientation(psi=0, theta=45, phi=0)
)

# Check shape type
print(matrix.shape_type)  # "sphere"
print(fiber.shape_type)   # "prolate_spheroid"

# Save to JSON
save_ellipsoids_json('composite.json', [matrix, fiber],
                     prop_names=['E', 'nu', 'alpha'])

# Load from JSON
phases = load_ellipsoids_json('composite.json')
```

### Layer Example (Laminates)

```python
from simcoon.solver.micromechanics import Layer, save_layers_json

# Create laminate layers
layer_0 = Layer(
    number=0,
    concentration=0.5,
    umat_name="ELISO",
    props=np.array([130000, 0.3, 1e-5]),
    geometry_orientation=GeometryOrientation(psi=0, theta=0, phi=0),  # 0-degree ply
    layerdown=1
)

layer_90 = Layer(
    number=1,
    concentration=0.5,
    umat_name="ELISO",
    props=np.array([130000, 0.3, 1e-5]),
    geometry_orientation=GeometryOrientation(psi=90, theta=0, phi=0),  # 90-degree ply
    layerup=0
)

save_layers_json('laminate.json', [layer_0, layer_90])
```

### JSON Format Examples

**ellipsoids.json**:
```json
{
  "ellipsoids": [
    {
      "number": 0,
      "coatingof": 0,
      "umat_name": "ELISO",
      "save": 1,
      "concentration": 0.7,
      "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
      "semi_axes": {"a1": 1, "a2": 1, "a3": 1},
      "geometry_orientation": {"psi": 0, "theta": 0, "phi": 0},
      "nstatev": 1,
      "props": {"E": 3500, "nu": 0.35, "alpha": 6e-5}
    }
  ]
}
```

**layers.json**:
```json
{
  "layers": [
    {
      "number": 0,
      "umat_name": "ELISO",
      "concentration": 0.5,
      "geometry_orientation": {"psi": 0, "theta": 90, "phi": -90},
      "props": {"E": 130000, "nu": 0.3, "alpha": 1e-5},
      "layerup": -1,
      "layerdown": 1
    }
  ]
}
```

---

## 5. API Comparison Tables

### Solver Functions

| v1.x | v2.0 | Notes |
|------|------|-------|
| `sim.solver(umat, props, ...)` | `Solver(blocks=[...]).solve()` | Class-based API |
| `sim.read_matprops(file)` | `load_material_json(file)` | JSON format |
| `sim.read_path(file)` | `load_path_json(file)` | Returns Block objects |
| File output parsing | `history[i].sigma`, etc. | Direct attribute access |

### Identification Functions

| v1.x | v2.0 | Notes |
|------|------|-------|
| `sim.identification()` | `IdentificationProblem` + optimizers | Full Python control |
| `sim.calc_cost()` | `mse()`, `mae()`, `r2()`, etc. | Multiple cost functions |
| File-based parameters | `ParameterSpec` class | Programmatic definition |

### Micromechanics I/O

| v1.x File | v2.0 Class | v2.0 JSON I/O |
|-----------|------------|---------------|
| `Nphases.dat` | `Phase` | `load_phases_json()` / `save_phases_json()` |
| `Nlayers.dat` | `Layer` | `load_layers_json()` / `save_layers_json()` |
| `Nellipsoids.dat` | `Ellipsoid` | `load_ellipsoids_json()` / `save_ellipsoids_json()` |
| `Ncylinders.dat` | `Cylinder` | `load_cylinders_json()` / `save_cylinders_json()` |

### Control Type Mappings

| v1.x Integer | v2.0 String | Description |
|--------------|-------------|-------------|
| 1 | `'small_strain'` | Infinitesimal strain |
| 2 | `'green_lagrange'` | Green-Lagrange strain |
| 3 | `'logarithmic'` | Logarithmic (Hencky) strain |
| 4 | `'biot'` | Biot strain |
| 5 | `'F'` | Deformation gradient control |
| 6 | `'gradU'` | Displacement gradient control |

### Corate Type Mappings

| v1.x Integer | v2.0 String | Description |
|--------------|-------------|-------------|
| 0 | `'jaumann'` | Jaumann (Zaremba-Jaumann) rate |
| 1 | `'green_naghdi'` | Green-Naghdi rate |
| 2 | `'logarithmic'` | Logarithmic rate |
| 3 | `'logarithmic_R'` | Logarithmic rate with rotation |
| 4 | `'truesdell'` | Truesdell rate |
| 5 | `'logarithmic_F'` | Logarithmic rate from F |

---

## 6. Common Migration Patterns

### Pattern 1: Simple Tension Test

**v1.x:**
```python
# Create path.txt manually, then:
sim.solver("ELISO", props, 1, 0, 0, 0, 0, 2, "data", "results", "path.txt", "output.txt")
data = np.loadtxt("results/output_global-0.txt")
```

**v2.0:**
```python
step = StepMeca(DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
                control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'])
block = Block(steps=[step], umat_name="ELISO", props=props, nstatev=1)
history = Solver(blocks=[block]).solve()
stress_strain = np.array([[h.Etot[0], h.sigma[0]] for h in history])
```

### Pattern 2: Cyclic Loading

**v2.0:**
```python
# Loading
step1 = StepMeca(
    DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress']
)
# Unloading
step2 = StepMeca(
    DEtot_end=np.array([-0.04, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress']
)
# Reloading
step3 = StepMeca(
    DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress']
)

block = Block(
    steps=[step1, step2, step3],
    umat_name="EPICP",
    props=props,
    nstatev=8,
    ncycle=5  # Repeat 5 times
)
```

### Pattern 3: Mixed Strain/Stress Control

**v2.0:**
```python
# Uniaxial stress with lateral strain measurement
step = StepMeca(
    DEtot_end=np.array([0, 0, 0, 0, 0, 0]),      # Strain targets (ignored for stress-controlled)
    Dsigma_end=np.array([500, 0, 0, 0, 0, 0]),   # Apply 500 MPa in direction 1
    control=['stress', 'stress', 'stress', 'stress', 'stress', 'stress']
)
```

### Pattern 4: Finite Strain with Deformation Gradient

**v2.0 (JSON path.json):**
```json
{
  "blocks": [{
    "type": "mechanical",
    "control_type": "deformation_gradient",
    "steps": [{
      "time": 5.0,
      "Dn_inc": 100,
      "F": [[5.0, 0, 0], [0, 0.447, 0], [0, 0, 0.447]]
    }]
  }]
}
```

### Pattern 5: Thermomechanical Loading

**v2.0:**
```python
from simcoon.solver import StepThermomeca

step = StepThermomeca(
    DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    DT_end=100.0,                    # Heat by 100 K
    thermal_control='temperature'    # Temperature-controlled (not adiabatic)
)

block = Block(
    steps=[step],
    umat_name="ELISO",
    umat_type="thermomechanical",
    props=props,
    nstatev=1
)
```

### Pattern 6: Convert Existing Text Files to JSON

```python
# Helper to migrate Nellipsoids.dat to JSON
def migrate_ellipsoids_file(txt_path, json_path):
    """Convert legacy Nellipsoids.dat to JSON format."""
    from simcoon.solver.micromechanics import Ellipsoid, save_ellipsoids_json

    ellipsoids = []
    with open(txt_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            # Parse according to Nellipsoids.dat format
            # (Adjust indices based on your actual file format)
            ell = Ellipsoid(
                number=int(parts[0]),
                coatingof=int(parts[1]),
                umat_name=parts[2],
                save=int(parts[3]),
                concentration=float(parts[4]),
                # ... parse remaining fields
            )
            ellipsoids.append(ell)

    save_ellipsoids_json(json_path, ellipsoids)
```

---

## 7. Troubleshooting

### Error: `AttributeError: module 'simcoon' has no attribute 'solver'`

**Cause:** Calling the removed `sim.solver()` function.

**Solution:** Use the new class-based API:
```python
from simcoon.solver import Solver, Block, StepMeca
# ... define steps and blocks ...
solver = Solver(blocks=[block])
history = solver.solve()
```

### Error: `TypeError: solver() missing required argument`

**Cause:** Mixing v1.x function signature with v2.0.

**Solution:** The new `Solver` class constructor only takes `blocks` and optional parameters:
```python
solver = Solver(blocks=[block], max_iter=10, tol=1e-9)
```

### Error: `ImportError: cannot import name 'identification' from 'simcoon'`

**Cause:** Using old import path.

**Solution:** Import from the new module:
```python
from simcoon.identification import IdentificationProblem, levenberg_marquardt
```

### Error: `ModuleNotFoundError: No module named 'scipy'`

**Cause:** scipy is required for the identification module.

**Solution:**
```bash
pip install scipy
```

### Error: `RuntimeError: Step failed to converge`

**Cause:** Newton-Raphson iteration did not converge within maximum sub-increments.

**Solution:**
1. Increase `Dn_inc` (maximum increments per step)
2. Decrease `Dn_init` (start with smaller increments)
3. Check material properties for physical consistency
4. Reduce strain/stress increment magnitude

```python
step = StepMeca(
    DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
    Dn_init=100,   # More increments
    Dn_inc=500,    # Allow more sub-increments
    Dn_mini=10     # Don't reduce below this
)
```

### Warning: Slow Solver Performance

**Cause:** History storage allocating too much memory.

**Solution:** The v2.0 solver uses lightweight `HistoryPoint` objects. For very long simulations, consider accessing results periodically:
```python
# Memory-efficient approach for very long simulations
solver = Solver(blocks=[block])
history = solver.solve()

# Process in chunks if needed
for i in range(0, len(history), 1000):
    chunk = history[i:i+1000]
    # Process chunk...
```

### Issue: JSON Files Not Loading Correctly

**Cause:** JSON syntax errors or missing required fields.

**Solution:** Validate JSON structure:
```python
import json

# Check JSON syntax
with open('material.json', 'r') as f:
    try:
        data = json.load(f)
        print("Valid JSON")
        print(json.dumps(data, indent=2))
    except json.JSONDecodeError as e:
        print(f"JSON Error: {e}")
```

Required fields for material.json:
- `name`: UMAT name (string)
- `props`: Material properties (object or array)
- `nstatev`: Number of state variables (integer)

### Issue: Results Don't Match v1.x

**Possible causes:**
1. Different control type/corate type mappings
2. Different increment sizes
3. Different convergence tolerances

**Solution:** Ensure equivalent settings:
```python
# v2.0 equivalent to v1.x solver_type=0, corate_type=2
block = Block(
    steps=[step],
    control_type='small_strain',  # was integer 1
    corate_type='logarithmic',    # was integer 2
    ...
)

solver = Solver(
    blocks=[block],
    tol=1e-9,         # Match v1.x tolerance
    max_iter=10       # Match v1.x iterations
)
```

---

## Additional Resources

- [Full Solver Documentation](simulation/solver.rst)
- [Micromechanics Documentation](simulation/micromechanics.rst)
- [Examples Repository](https://github.com/3MAH/simcoon/tree/master/examples)
- [Issue Tracker](https://github.com/3MAH/simcoon/issues)

---

*This migration guide covers Simcoon v2.0. For the latest updates, check the [official documentation](https://3mah.github.io/simcoon-docs/).*
