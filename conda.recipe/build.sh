#!/bin/bash

set -x

mkdir -p build
cd build
cmake .. -DCMAKE_INCLUDE_PATH=$PREFIX/include -DCMAKE_LIBRARY_PATH=$PREFIX/lib -DCMAKE_INSTALL_PREFIX=$PREFIX -Wno-dev
make
make install

cd $SRC_DIR/arma2numpy-builder
mkdir -p build
cd build
cmake .. -DCMAKE_INCLUDE_PATH=$PREFIX/include -DCMAKE_LIBRARY_PATH=$PREFIX/lib -DCMAKE_INSTALL_PREFIX=$PREFIX -Wno-dev -DCMAKE_BUILD_TYPE=Release
make
make install

cd $SRC_DIR/simcoon-python-builder
mkdir -p build
cd build
cmake .. -DCMAKE_INCLUDE_PATH=$PREFIX/include -DCMAKE_LIBRARY_PATH=$PREFIX/lib -DCMAKE_INSTALL_PREFIX=$PREFIX -Wno-dev -DCMAKE_BUILD_TYPE=Release
make
make install
cp lib/simmit.so $SRC_DIR/python-setup/simcoon/

cd $SRC_DIR/python-setup
pip install .

