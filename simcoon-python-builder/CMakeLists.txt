# CMakeLists.txt for Python bindings
# This is called from the parent CMakeLists.txt when building Python bindings
#
# Two modes:
#   1. SKBUILD (pip install): Bundle libsimcoon with simmit in wheel, RPATH=$ORIGIN
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

# Set RPATH based on build mode
if(CONDA_BUILD)
  # Conda: link to system libsimcoon in $PREFIX/lib
  if(APPLE)
    set_target_properties(simmit PROPERTIES
        INSTALL_RPATH "@loader_path/../../../"
        BUILD_WITH_INSTALL_RPATH ON
    )
  elseif(UNIX)
    set_target_properties(simmit PROPERTIES
        INSTALL_RPATH "$ORIGIN/../../../"
        BUILD_WITH_INSTALL_RPATH ON
    )
  endif()
else()
  # Wheel: libsimcoon is in same directory as simmit
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
endif()

# Set output file extension
if(MSVC)
  set_target_properties(simmit PROPERTIES PREFIX "" SUFFIX ".pyd")
else()
  set_target_properties(simmit PROPERTIES PREFIX "" SUFFIX ".so")
endif()

# Disable precompiled headers for this target (causes issues with carma includes)
set_target_properties(simmit PROPERTIES DISABLE_PRECOMPILE_HEADERS ON)

# Add include directories for Python wrapper headers
target_include_directories(simmit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link dependencies (simcoon provides Armadillo transitively)
target_link_libraries(simmit PRIVATE carma::carma Python3::Module simcoon)

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

  # Install simmit.so to site-packages/simcoon/
  install(TARGETS simmit
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
  install(TARGETS simmit
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
