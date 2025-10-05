# CMakeLists.txt for Python bindings
# This is called from the parent CMakeLists.txt when BUILD_PYTHON_BINDINGS is ON

cmake_minimum_required(VERSION 3.15)

# Option to control Python package installation location
option(SIMCOON_PYTHON_INSTALL_TO_SOURCE "Install Python package to source tree (development only)" OFF)

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

# Determine Python package installation directory
if(SIMCOON_PYTHON_INSTALL_TO_SOURCE)
  message(WARNING "Installing Python package to source tree - recommended for development only")
  set(PYTHON_INSTALL_DIR ${CMAKE_SOURCE_DIR}/python-setup/simcoon)
else()
  set(PYTHON_INSTALL_DIR ${CMAKE_BINARY_DIR}/python-package/simcoon)
  message(STATUS "Python package will be installed to build directory: ${PYTHON_INSTALL_DIR}")
endif()

# Install the Python extension module
install(TARGETS simmit
        DESTINATION ${PYTHON_INSTALL_DIR}
        COMPONENT python)

# Install the simcoon shared library to python-setup
install(FILES $<TARGET_FILE:simcoon>
        DESTINATION ${PYTHON_INSTALL_DIR}
        COMPONENT python)

# Install runtime dependency DLLs on Windows
# (DLL paths are found and stored by parent CMakeLists.txt)
if(WIN32 AND PYTHON_RUNTIME_DLLS)
  install(FILES ${PYTHON_RUNTIME_DLLS}
          DESTINATION ${PYTHON_INSTALL_DIR}
          COMPONENT python)
  list(LENGTH PYTHON_RUNTIME_DLLS DLL_COUNT)
  message(STATUS "Python package will include ${DLL_COUNT} runtime dependency DLLs")
endif()

################################################################################
# TESTS

add_test(NAME simmit_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/simmit_test
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/simmit_test/run_test.py)

# Set PYTHONPATH so the test can find the simcoon module
# Point to parent directory of PYTHON_INSTALL_DIR since it contains the 'simcoon' package directory
get_filename_component(PYTHON_PATH_DIR ${PYTHON_INSTALL_DIR} DIRECTORY)
if(WIN32)
  set(PATH_SEP ";")
else()
  set(PATH_SEP ":")
endif()
set_tests_properties(simmit_test PROPERTIES
  ENVIRONMENT "PYTHONPATH=${PYTHON_PATH_DIR}${PATH_SEP}$ENV{PYTHONPATH}"
)
