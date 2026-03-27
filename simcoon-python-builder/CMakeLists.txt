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

# All dependencies (Python, pybind11, carma, Armadillo, simcoon) are already found by parent

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
    # @loader_path          → find libsimcoon.dylib next to _core.so
    # @loader_path/../../.. → $CONDA_PREFIX/lib (conda's libomp.dylib)
    # /usr/local/lib        → CI / manual installs
    # /opt/homebrew/opt/libomp/lib → Homebrew on Apple Silicon
    set_target_properties(_core PROPERTIES
        INSTALL_RPATH "@loader_path;@loader_path/../../..;/usr/local/lib;/opt/homebrew/opt/libomp/lib"
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

# Use LLVM OpenMP on MSVC to avoid conflicts with the legacy vcomp140.dll
# runtime.  Requires CMake >= 3.30 and VS 2022 with the LLVM/Clang component.
# The compiler auto-links libomp via #pragma — no explicit path needed.
if(MSVC)
  set(OpenMP_RUNTIME_MSVC "llvm")
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(_core PRIVATE OpenMP::OpenMP_CXX)

  # Bundle the LLVM OpenMP DLL in the wheel on Windows
  if(WIN32 AND NOT CONDA_BUILD)
    # MSVC names the DLL by architecture: libomp140.<arch>.dll
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "ARM64|aarch64")
      set(_OMP_DLL_NAME "libomp140.aarch64.dll")
    else()
      set(_OMP_DLL_NAME "libomp140.x86_64.dll")
    endif()

    # Search in VS redistribution directories
    cmake_path(GET CMAKE_CXX_COMPILER PARENT_PATH _MSVC_BIN_DIR)
    find_file(_OMP_DLL ${_OMP_DLL_NAME}
      PATHS
        "${_MSVC_BIN_DIR}/../../Redist/MSVC" # VS Redist tree
        "ENV VCToolsRedistDir"
      PATH_SUFFIXES
        "debug_nonredist/x64/Microsoft.VC142.OpenMP.LLVM"
        "debug_nonredist/arm64/Microsoft.VC142.OpenMP.LLVM"
      NO_DEFAULT_PATH
    )

    if(_OMP_DLL)
      install(FILES "${_OMP_DLL}" DESTINATION . COMPONENT python)
      message(STATUS "Bundling OpenMP DLL: ${_OMP_DLL}")
    else()
      message(WARNING "Could not find ${_OMP_DLL_NAME} — wheel may require libomp at runtime")
    endif()
  endif()
endif()

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
