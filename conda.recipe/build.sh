#!/bin/bash

set -ex

cd $SRC_DIR

# Set pybind11 path for CMake
PYBIND11_DIR=$($PYTHON -c "import pybind11; print(pybind11.get_cmake_dir())")
if [ -z "$PYBIND11_DIR" ]; then
  echo "ERROR: Could not find pybind11 CMake directory"
  exit 1
fi
echo "Using pybind11 from: $PYBIND11_DIR"

# Configure with Python bindings enabled
cmake ${CMAKE_ARGS} -G Ninja -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX \
  -DCMAKE_INCLUDE_PATH=$PREFIX/include \
  -DCMAKE_LIBRARY_PATH=$PREFIX/lib \
  -DSIMCOON_BUILD_PYTHON_BINDINGS:BOOL=ON \
  -DSIMCOON_BUILD_TESTS:BOOL=ON \
  -DPython3_EXECUTABLE:FILEPATH=$PYTHON \
  -Dpybind11_DIR:PATH=$PYBIND11_DIR

# Build everything (C++ library and Python bindings)
cmake --build build -j${CPU_COUNT}

# Install C++ library to conda prefix
cmake --install build

# Install Python package components to build/python-package
cmake --install build --component python

# Install Python package from the generated package directory
$PYTHON -m pip install $SRC_DIR/build/python-package --no-deps -vv
