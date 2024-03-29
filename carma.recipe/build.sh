#!/bin/bash

set -ex

cd $SRC_DIR

cmake ${CMAKE_ARGS} -S . -B build \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX:path=$PREFIX \
  -D CARMA_INSTALL_LIB=ON

cmake --build build --config Release
cmake --install build
