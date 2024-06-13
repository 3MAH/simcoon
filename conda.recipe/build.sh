#!/bin/bash

set -ex

cd $SRC_DIR

cmake ${CMAKE_ARGS} -G Ninja -S . -B build \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX:path=$PREFIX \
  -D CMAKE_INCLUDE_PATH=$PREFIX/include \
  -D CMAKE_LIBRARY_PATH=$PREFIX/lib \

cmake -G Ninja --build build \ 
  -j${CPU_COUNT} --target simcoon --config Release
cmake -G Ninja --install build

cd $SRC_DIR/simcoon-python-builder

cmake ${CMAKE_ARGS} -G Ninja -S . -B build \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX:path=$PREFIX \
  -D CMAKE_INCLUDE_PATH=$PREFIX/include \
  -D CMAKE_LIBRARY_PATH=$PREFIX/lib

cmake -G Ninja --build build \
  -j${CPU_COUNT} --target simmit --config Release
cmake -G Ninja --install build
cp $SRC_DIR/simcoon-python-builder/build/lib/simmit.so $SRC_DIR/python-setup/simcoon/

cd $SRC_DIR/python-setup
$PYTHON -m pip install .
