# CMake generated Testfile for 
# Source directory: /scratch/source/simcoon/simcoon-python-builder
# Build directory: /scratch/source/simcoon/simcoon-python-builder/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Tarma2numpy "/usr/local/bin/python" "/scratch/source/simcoon/simcoon-python-builder/test/arma2numpy/run_test.py")
set_tests_properties(Tarma2numpy PROPERTIES  WORKING_DIRECTORY "/scratch/source/simcoon/simcoon-python-builder/test/arma2numpy")
