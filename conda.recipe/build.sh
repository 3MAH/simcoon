#!/bin/bash

set -ex

cd $SRC_DIR

# Set pybind11 path for CMake
PYBIND11_DIR=$($PYTHON -c "import pybind11; print(pybind11.get_cmake_dir())")
if [ -z "$PYBIND11_DIR" ]; then
  echo "ERROR: Could not find pybind11 CMake directory"
  exit 1
fi
echo "Using pybind11 from: $PYBIND11_DIR"

# The compiler activation may name the binutils without a path (gcc_linux-64 16.2.0 build 3);
# CMake then resolves them against the source directory. Pass them with their full path.
args=""
for arg in ${CMAKE_ARGS}; do
  case "$arg" in
    -DCMAKE_*AR=*|-DCMAKE_*RANLIB=*|-DCMAKE_LINKER=*|-DCMAKE_STRIP=*)
      tool="${arg#*=}"
      if [[ "$tool" != */* ]] && command -v "$tool" >/dev/null 2>&1; then
        arg="${arg%%=*}=$(command -v "$tool")"
      fi
      ;;
  esac
  args="$args $arg"
done
CMAKE_ARGS="$args"

# Configure with Python bindings enabled
# CONDA_BUILD env var is detected by CMakeLists.txt to enable SIMCOON_BUILD_PYTHON
cmake ${CMAKE_ARGS} -G Ninja -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX \
  -DCMAKE_INCLUDE_PATH=$PREFIX/include \
  -DCMAKE_LIBRARY_PATH=$PREFIX/lib \
  -DSIMCOON_BUILD_TESTS:BOOL=OFF \
  -DPython3_EXECUTABLE:FILEPATH=$PYTHON \
  -Dpybind11_DIR:PATH=$PYBIND11_DIR

# Build everything (C++ library and Python bindings)
cmake --build build -j${CPU_COUNT}

# Install C++ library to $PREFIX (lib, include, cmake)
# and Python bindings to site-packages
cmake --install build
