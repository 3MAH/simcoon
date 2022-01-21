#!/bin/bash

set -x

#Find the current directory
current_dir=$(pwd)

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
#cp build/lib/simmit.so $PREFIX/lib
cd ..

if [ $OS = "Mac" ]
then
    install_name_tool -change libsimcoon.dylib @rpath/libsimcoon.dylib $PREFIX/lib/libarma2numpy.dylib
    install_name_tool -change ${current_dir}/simcoon-python-builder/build/lib/libarma2numpy.dylib  @rpath/libarma2numpy.dylib $PREFIX/lib/python${python_version}/site-packages/simcoon/simmit.so
    install_name_tool -change libsimcoon.dylib @rpath/libsimcoon.dylib $PREFIX/lib/python${python_version}/site-packages/simcoon/simmit.so
fi

cd python-setup
pip install .

