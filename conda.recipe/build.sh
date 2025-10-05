#!/bin/bash

set -ex

cd $SRC_DIR

# Set pybind11 path for CMake
PYBIND11_DIR=$($PYTHON -c "import pybind11; print(pybind11.get_cmake_dir())")

# Configure with Python bindings enabled
cmake ${CMAKE_ARGS} -G Ninja -S . -B build \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
  -D CMAKE_INCLUDE_PATH=$PREFIX/include \
  -D CMAKE_LIBRARY_PATH=$PREFIX/lib \
  -D SIMCOON_BUILD_PYTHON_BINDINGS:BOOL=ON \
  -D SIMCOON_BUILD_TESTS:BOOL=ON \
  -D Python3_EXECUTABLE:FILEPATH=$PYTHON \
  -D pybind11_DIR:PATH=$PYBIND11_DIR \
  -D CMAKE_PREFIX_PATH=$PYBIND11_DIR

# Build everything (C++ library and Python bindings)
cmake --build build -j${CPU_COUNT}

# Install C++ library
cmake --install build

# Install Python package components
cmake --install build --component python

# Install Python package
cd $SRC_DIR/build/python-package
$PYTHON -m pip install . --no-deps -vv
