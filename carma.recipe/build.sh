#!/bin/bash

set -x

mkdir -p build
cd build
cmake -DCARMA_INSTALL_LIB=ON -DCMAKE_INSTALL_PREFIX=$PREFIX ..
cmake --build . --config Release --target install
