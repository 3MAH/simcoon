#!/bin/bash

set -ex

cd $SRC_DIR

cmake ${CMAKE_ARGS} -S . -B build \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX:path=$PREFIX
  -D CMAKE_INCLUDE_PATH=$PREFIX/include
  -D CMAKE_LIBRARY_PATH=$PREFIX/lib

cmake --build build --config Release
cmake --install build

cd $SRC_DIR/simcoon-python-builder

cmake ${CMAKE_ARGS} -S . -B build \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX:path=$PREFIX
  -D CMAKE_INCLUDE_PATH=$PREFIX/include
  -D CMAKE_LIBRARY_PATH=$PREFIX/lib

cmake --build build --config Release
cmake --install build

mkdir -p build
cp lib/simmit.so $SRC_DIR/python-setup/simcoon/

cd $SRC_DIR/python-setup
pip install .
