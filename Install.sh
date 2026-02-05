#!/bin/bash

display_help() {
   echo ""
   echo "Usage: $0 -t test -n ncpus"
   echo "\t-t run the script without executing tests"
   echo "\t-n Number of cpus for the compilation"
   exit 1 # Exit script after printing help
}

test=1
ncpus=4

while getopts "tn:" opt
do
   case $opt in
      t ) test=0 ;;
      n ) ncpus="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$test" ] || [ -z $ncpus ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct

echo "Start of Simcoon compilation."

cp external/umat_plugin_ext.cpp testBin/Umats/UMEXT/external/umat_plugin_ext.cpp
cp external/umat_plugin_aba.cpp testBin/Umats/UMABA/external/umat_plugin_aba.cpp

blue=`tput setaf 4`
red=`tput setaf 1`
reset=`tput sgr0`

#Find the current directory
current_dir=$(pwd)

# Test if build exists and create it if necessary
if [ ! -d "${current_dir}/build" ]; then
    mkdir -p "${current_dir}/build"
    echo "Build folder created."
else
    echo "Build directory already exists."

    while true; do
        read -p "Do you want to erase old compilation files (Recommended: No)? " yn
        case $yn in
            [YyOo]* )
                # Remove all regular files, directories, and shared libraries recursively that are writable by you
                find "${current_dir}/build" -type f -name "*.so" -user "$USER" -exec rm -f {} \;
                find "${current_dir}/build" -mindepth 1 -user "$USER" -exec rm -rf {} \;
                echo "Old compilation files erased."
                break
                ;;
            [Nn]* )
                echo "Keeping existing build files."
                break
                ;;
            * )
                echo "Please answer yes (y) or no (n)."
                ;;
        esac
    done
fi

#Build Simcoon
uv run cmake -B build -G Ninja \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_PREFIX_PATH="$(uv run python -m pybind11 --cmakedir)"
cmake --build build -j${ncpus}
cmake --install build --prefix "${CONDA_PREFIX}"

uv pip install ./build/python-package
uv run python -c "import simcoon; import simcoon.simmit; print('simcoon version:', simcoon.__version__)"

if [ $test -eq 1 ]
then
    uv run ctest -j${ncpus} --test-dir build --output-on-failure
fi

