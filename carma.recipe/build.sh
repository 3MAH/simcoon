#!/bin/bash

set -x

mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:path=$CONDA_PREFIX -DCARMA_INSTALL_LIB=ON
make
make install
