#!/bin/bash

set -x

mkdir -p build
cd build
cmake .. -DCMAKE_INCLUDE_PATH=$PREFIX/include -DCMAKE_LIBRARY_PATH=$PREFIX/lib -DCMAKE_INSTALL_PREFIX=$PREFIX -Wno-dev
make
make install
cd ..

cd simcoon-python-builder
mkdir -p build
cd build
cmake .. -DCMAKE_INCLUDE_PATH=$PREFIX/include -DCMAKE_LIBRARY_PATH=$PREFIX/lib -DCMAKE_INSTALL_PREFIX=$PREFIX -Wno-dev -DCMAKE_BUILD_TYPE=Release
make
make install
cd ..
cp build/lib/simmit.so ../python-setup/simcoon/
cp build/lib/simmit.so $PREFIX/lib
cd ..

cd python-setup
pip install .
