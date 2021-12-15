#!/bin/bash

# conda install -c conda-forge armadillo boost cgal numpy -y

anacondaloc=$(dirname "$(dirname "$(which python)")")

mkdir -p build
cd build
cmake .. -DCMAKE_INCLUDE_PATH=${anacondaloc}/include -DCMAKE_LIBRARY_PATH=${anacondaloc}/lib -DCMAKE_INSTALL_PREFIX=${anacondaloc} -Wno-dev
make
make install
cd ..

cd simcoon-python-builder
mkdir -p build
cd build
cmake .. -DCMAKE_INCLUDE_PATH=${anacondaloc}/include -DCMAKE_LIBRARY_PATH=${anacondaloc}/lib -Wno-dev -DCMAKE_BUILD_TYPE=Release
make
cd ..
cp -r include/* ${anacondaloc}/include
cp build/lib/libarma2numpy.so ${anacondaloc}/lib
cp build/lib/simmit.so ../python-setup/simcoon/
cd ..

cd python-setup/
${PYTHON} setup.py install