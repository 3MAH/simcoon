Simcoon
=========

[![Simcoon Logo](https://github.com/3MAH/3mah.github.io/blob/master/assets/images/logo_simcoon/simcoon_logo_text.png?raw=true)](https://github.com/3MAH/simcoon)

About
-----

Simcoon is a free, open-source library for the simulation of multiphysics systems. Its primarily objective was the developement of constitutive models for the simulation of heterogeneous materials, but now goes beyond with tools to facilitate their full-field simulation. Together with [microgen](https://github.com/3MAH/microgen) for the CAD and meshing of heterogeneous materials and [fedoo](https://github.com/3MAH/fedoo) our Finite Element solver, we offer a comprehensive simulation set for the in-depth analysis of heterogeneous materials.

Simcoon is developed with the aim to be a high-quality scientific library to facilitate the analysis of the complex, non-linear behavior of systems. It integrates tools to simulate the response of composite material response and thus integrates several algorithms for the analysis of heterogeneous materials.

Simcoon integrates
- a easy way to handle geometrical non-linearities : Use of Lagrangian measures, Eulerian measures and cumulative strains considering several spins : Jaumann, Green-Naghdi, Xi-Meyers-Bruhns logarithmic. With this last measure, cumulative strain correspond to a logarithmic strain measure and is the standard measure utilized for our constitutive laws.

Simcoon is a C++ library with emphasis on speed and ease-of-use, that offers a python interface to facilitate its use. Its principle focus is to provide tools to facilitate the implementation of up-to-date constitutive model for materials in Finite Element Analysis Packages. This is done by providing a C++ API to generate user material subroutine based on a library of functions. Also, Simconnn provides tools to analyse the behavior of material, considering loading at the material point level. Such tools include a thermomechanical solver, a software to predict effective properties of composites, and a built-in identification software (using a combined genetic-gradient based algorithm)

Simcoon is mainly developed by faculty and researchers from University of Bordeaux and the I2M Laboratory (Institut de d'Ingénierie et de Mécanique). Fruitful contribution came from the LEM3 laboratory in Metz, France, TU Bergakademie Freiberg in Germany and the TIMC-IMAG laboratory in Grenoble, France. It is released under the GNU General Public License: GPL, version 3.

[![GitHub license](https://img.shields.io/badge/licence-GPL%203-blue.svg)](https://github.com/chemiskyy/simcoon/blob/master/LICENSE.txt)

Simcoon make use and therefore include the FTensor library (http://www.wlandry.net/Projects/FTensor) for convenience. FTensor is a library that handle complex tensor computations. FTensor is released under the GNU General Public License: GPL, version 2. You can get it there (but is is already included in simcoon): (https://bitbucket.org/wlandry/ftensor)

Quick Start (Python)
--------------------

Simcoon v2.0 introduces a new Python-based solver API that replaces the legacy file-based C++ solver:

```python
import numpy as np
from simcoon.solver import Solver, Block, StepMeca

# Material properties: E, nu, alpha
props = np.array([210000.0, 0.3, 1e-5])

# Define uniaxial tension loading
step = StepMeca(
    DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),  # 1% strain
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=50
)

# Create simulation block
block = Block(
    steps=[step],
    umat_name='ELISO',
    props=props,
    nstatev=1
)

# Run simulation
solver = Solver(blocks=[block])
history = solver.solve()

# Extract results
strain = np.array([h.Etot[0] for h in history])
stress = np.array([h.sigma[0] for h in history])
```

For more examples, see `examples/umats/` and the [documentation](https://3mah.github.io/simcoon-docs/).

Documentation
--------------

Provider      | Status
--------      | ------
Documentation | [![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-blue?style=plastic)](https://3mah.github.io/simcoon-docs/)


Building doc :
requires doxygen, sphinx, breathe

```bash
conda install -c conda-forge doxygen -y && pip install sphinx sphinx-rtd-theme breathe
```

```bash
cd doxdocs && make html
```

open _build/index.html


Installation
------------

### Option 1: Install from Conda (Recommended)

The simplest way to install simcoon is directly with conda:

```bash
conda install -c conda-forge -c set3mah simcoon
```

In case of conflicts, create a new conda environment:
```bash
conda create --name simcoon_env
conda activate simcoon_env
conda install -c conda-forge -c set3mah simcoon
```

### Option 2: Build from Source

#### Prerequisites

Create and activate a conda environment:
```bash
conda create --name simcoon_build
conda activate simcoon_build
```

Install required dependencies:
```bash
# Compilers and build tools
conda install -c conda-forge cxx-compiler fortran-compiler cmake ninja

# Libraries
conda install -c conda-forge armadillo pybind11 numpy gtest carma

# Python testing
pip install pytest
```

For x86 architectures, you may also need MKL:
```bash
conda install -c conda-forge mkl
```

#### Build Instructions (without conda)

1. Clone or download the repository:
```bash
git clone https://github.com/3MAH/simcoon.git
cd simcoon
```

2. Install required dependencies using your system's package manager.

- On Debian/Ubuntu:

```bash
sudo apt-get install libarmadillo-dev libgtest-dev ninja-build
```

- On macOS with Homebrew:

```bash
brew install armadillo googletest
```

- On Windows with vcpkg:

```powershell
vcpkg install armadillo gtest
```

3. Configure and build the project:

**For Python users** (recommended):
```bash
pip install .
```

**For C++ development:**
```bash
# Configure and build
cmake -S . -B build -G Ninja -D CMAKE_BUILD_TYPE=Release
cmake --build build

# Run C++ tests
ctest --test-dir build --output-on-failure
```

#### Development Workflow

For active development with both C++ and Python:

```bash
# Install build dependencies first
uv pip install scikit-build-core pybind11 numpy

# Editable install (uv applies --no-build-isolation automatically via pyproject.toml)
uv pip install -e .[dev]

# After modifying C++ files, rebuild directly
cmake --build build/cp*

# Python changes take effect immediately (no rebuild needed)
```

The editable install creates a build directory at `build/{wheel_tag}` (e.g., `build/cp312-cp312-linux_x86_64`). The `[tool.uv]` config in `pyproject.toml` disables build isolation for simcoon, ensuring the CMake cache references your actual Python environment, enabling direct `cmake --build` commands for incremental rebuilds.

**Auto-rebuild on import**: Importing simcoon will automatically trigger a cmake rebuild if C++ files have changed:
```bash
python -c "import simcoon"  # Rebuilds if needed
uv run python -c "import simcoon"  # Also works
```

**Note**: If you add new C++ source files, re-run `uv pip install -e .[dev]` to reconfigure.

#### Build Options

- `SIMCOON_BUILD_TESTS` (default: ON) - Build C++ tests (CMake only)

#### Notes for macOS

For numpy versions earlier than 1.26.4 using the Accelerate framework:
```bash
pip install cython
pip install --no-binary :all: numpy
```


Authors
-------
* [Yves Chemisky](https://github.com/chemiskyy)
