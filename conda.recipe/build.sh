#!/bin/bash

set -ex

cd $SRC_DIR

# Configure C++ library and tests (Python bindings built separately via pip)
cmake ${CMAKE_ARGS} -G Ninja -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX \
  -DCMAKE_INCLUDE_PATH=$PREFIX/include \
  -DCMAKE_LIBRARY_PATH=$PREFIX/lib \
  -DSIMCOON_BUILD_PYTHON_BINDINGS:BOOL=OFF \
  -DSIMCOON_BUILD_TESTS:BOOL=ON

# Build C++ library
cmake --build build -j${CPU_COUNT}

# Install C++ library and headers to conda prefix
cmake --install build

# Install Python package via scikit-build-core
cd $SRC_DIR
$PYTHON -m pip install . --no-deps --no-build-isolation -vv
