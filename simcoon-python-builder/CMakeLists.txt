# CMakeLists.txt for Python bindings
# This is called from the parent CMakeLists.txt when BUILD_PYTHON_BINDINGS is ON

cmake_minimum_required(VERSION 3.15)

# Option to control Python package installation location
option(SIMCOON_PYTHON_INSTALL_TO_SOURCE "Install Python package to source tree (development only)" OFF)

# All dependencies (Python, pybind11, carma, Armadillo, simcoon) are already found by parent

################################################################################
# SIMCOON PYBIND11 MODULE (simmit)

pybind11_add_module(simmit
  # Main wrapper
  src/python_wrappers/boostpython_smartplus_wrappers.cpp

  # Continuum Mechanics
  src/python_wrappers/Libraries/Continuum_mechanics/constitutive.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/contimech.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/criteria.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/damage.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/hyperelastic.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/kinematics.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/Leff.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/objective_rates.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/recovery_props.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/stress.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/transfer.cpp
  src/python_wrappers/Libraries/Continuum_mechanics/umat.cpp

  # Homogenization
  src/python_wrappers/Libraries/Homogenization/eshelby.cpp

  # Identification
  src/python_wrappers/Libraries/Identification/constants.cpp
  src/python_wrappers/Libraries/Identification/identification.cpp
  src/python_wrappers/Libraries/Identification/optimize.cpp
  src/python_wrappers/Libraries/Identification/parameters.cpp

  # Material
  src/python_wrappers/Libraries/Material/ODF.cpp

  # Maths
  src/python_wrappers/Libraries/Maths/lagrange.cpp
  src/python_wrappers/Libraries/Maths/rotation.cpp

  # Solver
  src/python_wrappers/Libraries/Solver/read.cpp
  src/python_wrappers/Libraries/Solver/solver.cpp
)

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

# C++ standard inherited from parent (C++20)

# Add include directories for Python wrapper headers
target_include_directories(simmit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link dependencies (simcoon provides Armadillo transitively)
target_link_libraries(simmit PRIVATE carma::carma Python3::Module simcoon)

# Determine Python package installation directory
if(SIMCOON_PYTHON_INSTALL_TO_SOURCE)
  message(WARNING "Installing Python package to source tree - recommended for development only")
  set(PYTHON_PACKAGE_ROOT ${CMAKE_SOURCE_DIR}/python-setup)
  set(PYTHON_INSTALL_DIR ${CMAKE_SOURCE_DIR}/python-setup/simcoon)
else()
  set(PYTHON_PACKAGE_ROOT ${CMAKE_BINARY_DIR}/python-package)
  set(PYTHON_INSTALL_DIR ${CMAKE_BINARY_DIR}/python-package/simcoon)
  message(STATUS "Python package will be installed to build directory: ${PYTHON_PACKAGE_ROOT}")
endif()

# Install setup.py and Python source files to create a complete package
install(FILES ${CMAKE_SOURCE_DIR}/python-setup/setup.py
        DESTINATION ${PYTHON_PACKAGE_ROOT}
        COMPONENT python)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/python-setup/simcoon/
        DESTINATION ${PYTHON_INSTALL_DIR}
        COMPONENT python
        FILES_MATCHING PATTERN "*.py")

# Install the Python extension module
install(TARGETS simmit
        DESTINATION ${PYTHON_INSTALL_DIR}
        COMPONENT python)

# Install the simcoon shared library
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
