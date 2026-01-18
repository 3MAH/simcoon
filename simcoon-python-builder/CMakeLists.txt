# CMakeLists.txt for Python bindings
# This is called from the parent CMakeLists.txt when SIMCOON_BUILD_PYTHON_BINDINGS is ON

# Option to control Python package installation location
option(SIMCOON_PYTHON_INSTALL_TO_SOURCE "Install Python package to source tree (development only)" OFF)
project(simcoon)

#set(CMAKE_VERBOSE_MAKEFILE ON)

# The version number
set(VERSION_MAJOR 1)
set(VERSION_MINOR 9)
set(VERSION_PATCH 4)

if(VERSION_PATCH MATCHES "0")
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}")
else()
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()

set(CMAKE_DISABLE_PRECOMPILE_HEADERS ON)

##### RPATH HANDLING ####
#set(CMAKE_MACOSX_RPATH 1)
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
################


# put our local cmake find scripts at the beginning of the cmake
# module search path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (MSVC) #if windows using visual c++

  # Specify the output directory for each configuration (Debug, Release, etc.)
  foreach(CONFIG_TYPE ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG_TYPE} CONFIG)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/bin")
  endforeach()
endif()

#Set executable files and library files
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

find_package(Python3 COMPONENTS Interpreter Development NumPy)

#Inclusion of Armadillo
if (MSVC)
  find_package(pybind11 REQUIRED)
  find_package(carma CONFIG REQUIRED)
  find_package(blas REQUIRED)
  find_package(lapack REQUIRED)
  set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})  
else()
  find_package(pybind11 REQUIRED)
  find_package(Armadillo 12.6 REQUIRED)
  find_package(carma CONFIG REQUIRED)
endif()

if (NOT DEFINED USE_OPENMP)
  set(USE_OPENMP True)
endif()

if (USE_OPENMP)
  find_package(OpenMP)
endif()
find_package(simcoon REQUIRED)

#${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}

include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})
include_directories(SYSTEM ${SIMCOON_INCLUDE_DIRS})
include_directories(SYSTEM ${CARMA_INCLUDE_DIR})

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "SIMCOON_INCLUDE_DIRS = ${SIMCOON_INCLUDE_DIRS}")
message(STATUS "SIMCOON_LIBRARIES = ${SIMCOON_LIBRARIES}")
message(STATUS "ARMADILLO_INCLUDE_DIRS = ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "ARMADILLO_LIBRARIES = ${ARMADILLO_LIBRARIES}")
message(STATUS "carma_INCLUDE_DIR = ${carma_INCLUDE_DIR}")

# Build type
if(NOT CMAKE_BUILD_TYPE)  # Debug by default
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel"
        FORCE)
endif()

# Set the CFLAGS and CXXFLAGS depending on the options the user specified.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-unused-parameter")


if (MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
  if (USE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict= /openmp:llvm /Y-") #\Y- to disable precompile header (don't work for an unknown reason)
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict= /openmp /Y-") #\Y- to disable precompile header (don't work for an unknown reason)
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict=")
  endif()
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
elseif (UNIX AND NOT APPLE)
  if (DEBUG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O0 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++20")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20")
  else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O3 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++20")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20")
  endif()
else()
  if (DEBUG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -stdlib=libc++ -O0 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
  else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -stdlib=libc++ -O3 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
  endif()
endif()

#Inclusion of public headers
include_directories(include)

#Command file to get all the files in the src/ and include/
file(GLOB_RECURSE source_files_simmit src/python_wrappers/* include/simcoon/python_wrappers/*)
file(GLOB_RECURSE test_files_carma test/Carma/*.cpp test/Carma/*.hpp)
file(GLOB_RECURSE test_files_exception test/exception_test/*.cpp test/exception_test/*.hpp)

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

# Disable precompiled headers for this target (causes issues with carma includes)
set_target_properties(simmit PROPERTIES DISABLE_PRECOMPILE_HEADERS ON)

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
if (MSVC) #if windows using visual c++

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/exception_test)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/exception_test)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/exception_test)

  # Specify the output directory for each configuration (Debug, Release, etc.)
  foreach(CONFIG_TYPE ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG_TYPE} CONFIG)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_CURRENT_SOURCE_DIR}/test/exception_test")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_CURRENT_SOURCE_DIR}/test/exception_test")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_CURRENT_SOURCE_DIR}/test/exception_test")
  endforeach()
endif()

set(LIBRARY_OUTPUT_PATH {CMAKE_CURRENT_SOURCE_DIR}/test)
pybind11_add_module(test_exception ${test_files_exception})

# Set output directories based on the build type
if(MSVC)
  set_target_properties(test_exception PROPERTIES
    PREFIX ""
    SUFFIX ".pyd"
  )
else()
  # For non-MSVC platforms
  set_target_properties(test_exception PROPERTIES
    PREFIX ""
    SUFFIX ".so"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test/exception_test"
  )
endif()
#Link libraries
target_link_libraries(test_exception PRIVATE carma::carma ${ARMADILLO_LIBRARIES} Python3::Module)

#Finally add it to test execution -
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME test_simcoon_exception WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/exception_test
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/exception_test/run_test_exception.py)

#Finally add it to test execution -
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME simmit_test WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/simmit_test
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/simmit_test/run_test.py)

# Install the Python extension module
install(TARGETS simmit
        DESTINATION ${PYTHON_INSTALL_DIR}
        COMPONENT python)

# Install the simcoon shared library with SOVERSION symlinks
# Must use install(TARGETS) to get proper symlink handling
install(TARGETS simcoon
        LIBRARY DESTINATION ${PYTHON_INSTALL_DIR}
          COMPONENT python
          NAMELINK_COMPONENT python
        RUNTIME DESTINATION ${PYTHON_INSTALL_DIR}
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
