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

Documentation
--------------

Provider      | Status
--------      | ------
Github | [![Documentation Status](https://readthedocs.org/projects/simcoon/badge/?version=latest&style=plastic)](https://3mah.github.io/simcoon-docs/)


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

For Python 3.11 or later, you may also need to install libboost:
```bash
conda install -c conda-forge libboost
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
conda install -c conda-forge armadillo boost pybind11 numpy gtest carma

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
sudo apt-get install libarmadillo-dev libboost-all-dev libgtest-dev ninja-build
```

- On macOS with Homebrew:

```bash
brew install armadillo boost googletest
```

- On Windows with vcpkg:

```powershell
vcpkg install armadillo boost-config boost-dll gtest
```

3. Configure and build the project:

**Linux/macOS:**
```bash
# Configure
cmake -S . -B build -G Ninja -D CMAKE_BUILD_TYPE=Release

# Build
cmake --build build

# Install Python package
pip install ./build/python-package
```

**Windows:**
```powershell
# Configure
cmake -S . -B build

# Build
cmake --build build --config Release

# Install Python package
pip install ./build/python-package
```

3. (Optional) Run tests:
```bash
ctest --test-dir build --output-on-failure
```

#### Build Options

- `SIMCOON_BUILD_PYTHON_BINDINGS` (default: ON) - Build Python bindings
- `SIMCOON_BUILD_TESTS` (default: ON) - Build C++ tests

To build only the C++ library without Python bindings:
```bash
cmake -S . -B build -D SIMCOON_BUILD_PYTHON_BINDINGS=OFF
cmake --build build
cmake --install build --prefix /path/to/install
```

#### Notes for macOS

For numpy versions earlier than 1.26.4 using the Accelerate framework:
```bash
pip install cython
pip install --no-binary :all: numpy
```


Authors
-------
* [Yves Chemisky](https://github.com/chemiskyy)
