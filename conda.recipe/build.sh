#!/bin/bash

set -ex

cd $SRC_DIR

# Configure with Python bindings enabled
cmake ${CMAKE_ARGS} -G Ninja -S . -B build \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
  -D CMAKE_INCLUDE_PATH=$PREFIX/include \
  -D CMAKE_LIBRARY_PATH=$PREFIX/lib

# Build everything (C++ library and Python bindings)
cmake --build build -j${CPU_COUNT}

# Install C++ library
cmake --install build

# Install Python package components
cmake --install build --component python

# Install Python package
cd $SRC_DIR/build/python-package
$PYTHON -m pip install . --no-deps -vv
