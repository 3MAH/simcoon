# CMakeLists.txt for Python bindings
# This is called from the parent CMakeLists.txt when building Python bindings
#
# Two modes:
#   1. SKBUILD (pip install): Bundle libsimcoon with _core in wheel, RPATH=$ORIGIN
#   2. CONDA_BUILD (conda-build): Link to system libsimcoon, install to site-packages

if(DEFINED SKBUILD)
  message(STATUS "Building Python bindings for wheel (scikit-build-core)")
elseif(CONDA_BUILD)
  message(STATUS "Building Python bindings for conda")
else()
  message(STATUS "Building Python bindings (standalone)")
endif()

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
# SIMCOON PYBIND11 MODULE (_core)

pybind11_add_module(_core
  # Main wrapper
  src/python_wrappers/python_module.cpp

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

# Set RPATH based on build mode
if(CONDA_BUILD)
  # Conda: link to system libsimcoon in $PREFIX/lib
  if(APPLE)
    set_target_properties(_core PROPERTIES
        INSTALL_RPATH "@loader_path/../../../"
        BUILD_WITH_INSTALL_RPATH ON
    )
  elseif(UNIX)
    set_target_properties(_core PROPERTIES
        INSTALL_RPATH "$ORIGIN/../../../"
        BUILD_WITH_INSTALL_RPATH ON
    )
  endif()
else()
  # Wheel: libsimcoon is in same directory as _core
  if(APPLE)
    set_target_properties(_core PROPERTIES
        INSTALL_RPATH "@loader_path"
        BUILD_WITH_INSTALL_RPATH ON
    )
  elseif(UNIX)
    set_target_properties(_core PROPERTIES
        INSTALL_RPATH "$ORIGIN"
        BUILD_WITH_INSTALL_RPATH ON
    )
  endif()
endif()

# Set output file extension
if(MSVC)
  set_target_properties(_core PROPERTIES PREFIX "" SUFFIX ".pyd")
else()
  set_target_properties(_core PROPERTIES PREFIX "" SUFFIX ".so")
endif()

# Disable precompiled headers for this target (causes issues with carma includes)
set_target_properties(_core PROPERTIES DISABLE_PRECOMPILE_HEADERS ON)

# Add include directories for Python wrapper headers
target_include_directories(_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link dependencies (simcoon provides Armadillo transitively)
target_link_libraries(_core PRIVATE carma::carma Python3::Module simcoon)

################################################################################
# INSTALL

if(CONDA_BUILD)
  # Conda: Install to site-packages/simcoon/
  # Get Python site-packages path
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('purelib'))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES_RAW
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  # Convert Windows backslashes to forward slashes to avoid CMake escape issues
  file(TO_CMAKE_PATH "${PYTHON_SITE_PACKAGES_RAW}" PYTHON_SITE_PACKAGES)
  message(STATUS "Python site-packages: ${PYTHON_SITE_PACKAGES}")

  # Install _core.so to site-packages/simcoon/
  install(TARGETS _core
          LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}/simcoon
          RUNTIME DESTINATION ${PYTHON_SITE_PACKAGES}/simcoon)

  # Install pure Python files to site-packages/simcoon/
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/python-setup/simcoon/
          DESTINATION ${PYTHON_SITE_PACKAGES}/simcoon
          FILES_MATCHING PATTERN "*.py")

  # Note: libsimcoon is already installed to $PREFIX/lib by parent CMakeLists.txt

else()
  # Wheel: scikit-build-core handles the wheel packaging

  # Install Python extension module to wheel
  install(TARGETS _core
          DESTINATION .
          COMPONENT python)

  # Install libsimcoon shared library to wheel
  install(TARGETS simcoon
          LIBRARY DESTINATION .
            COMPONENT python
            NAMELINK_COMPONENT python
          RUNTIME DESTINATION .
            COMPONENT python)

  # Install runtime dependency DLLs on Windows
  if(WIN32 AND PYTHON_RUNTIME_DLLS)
    install(FILES ${PYTHON_RUNTIME_DLLS}
            DESTINATION .
            COMPONENT python)
    list(LENGTH PYTHON_RUNTIME_DLLS DLL_COUNT)
    message(STATUS "Python package will include ${DLL_COUNT} runtime dependency DLLs")
  endif()
endif()

# Python tests are run separately via pytest
# See: pytest simcoon-python-builder/test/ -v
