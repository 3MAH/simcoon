# CMakeLists.txt for Python bindings
# This is called from the parent CMakeLists.txt when BUILD_PYTHON is ON

cmake_minimum_required(VERSION 3.15)

# All dependencies (Python, pybind11, carma, Armadillo, simcoon) are already found by parent

# Add this subproject's include directory
include_directories(include)

# Get source files for Python bindings
file(GLOB_RECURSE source_files_simmit src/python_wrappers/* include/simcoon/python_wrappers/*)

################################################################################
# SIMCOON PYBIND11 MODULE (simmit)

pybind11_add_module(simmit ${source_files_simmit})

# Set RPATH so simmit.so can find libsimcoon.so in the same directory
if(APPLE)
  set_target_properties(simmit PROPERTIES
      INSTALL_RPATH "@loader_path"
      BUILD_WITH_INSTALL_RPATH ON
  )
elseif(UNIX)
  set_target_properties(simmit PROPERTIES
      INSTALL_RPATH "$ORIGIN"
      BUILD_WITH_INSTALL_RPATH ON
  )
endif()

# Set output file extension
if(MSVC)
  set_target_properties(simmit PROPERTIES PREFIX "" SUFFIX ".pyd")
else()
  set_target_properties(simmit PROPERTIES PREFIX "" SUFFIX ".so")
endif()
set_target_properties(simmit PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

target_link_libraries(simmit PRIVATE ${ARMADILLO_LIBRARIES} carma::carma Python3::Module simcoon)

# Install the Python extension module
install(TARGETS simmit
        DESTINATION ${CMAKE_SOURCE_DIR}/python-setup/simcoon
        COMPONENT python)

# Install the simcoon shared library to python-setup
install(FILES $<TARGET_FILE:simcoon>
        DESTINATION ${CMAKE_SOURCE_DIR}/python-setup/simcoon
        COMPONENT python)

# Note: Dependency DLLs are installed by the main CMakeLists.txt

################################################################################
# TESTS

add_test(NAME simmit_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/simmit_test
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/simmit_test/run_test.py)

# Set PYTHONPATH so the test can find the simcoon module
set_tests_properties(simmit_test PROPERTIES
  ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/python-setup:$ENV{PYTHONPATH}"
)
