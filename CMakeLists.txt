# CMakeLists.txt

cmake_minimum_required(VERSION 3.15)

# Read version from pyproject.toml ([project].version)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml" PYPROJECT_FILE)
string(REGEX MATCH "\\[project\\][^[]*" PYPROJECT_PROJECT_SECTION "${PYPROJECT_FILE}")
string(REGEX MATCH "\nversion[ \t]*=[ \t]*\"([^\"]+)\"" _ "${PYPROJECT_PROJECT_SECTION}")
if(NOT CMAKE_MATCH_1)
  message(FATAL_ERROR "Could not extract [project].version from pyproject.toml")
endif()
set(SIMCOON_VERSION "${CMAKE_MATCH_1}")

# CMake project() version must be numeric; derive ABI version from PEP 440 version.
string(REGEX MATCH "^([0-9]+\\.[0-9]+\\.[0-9]+)" _ "${SIMCOON_VERSION}")
if(NOT CMAKE_MATCH_1)
  message(FATAL_ERROR "Could not derive numeric CMake version (X.Y.Z) from '${SIMCOON_VERSION}'")
endif()
set(SIMCOON_ABI_VERSION "${CMAKE_MATCH_1}")

project(simcoon VERSION ${SIMCOON_ABI_VERSION} LANGUAGES C CXX)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)

# Print the project version
message(STATUS "Project version: ${SIMCOON_VERSION} (ABI: ${PROJECT_VERSION})")

# Optionally use this version in packaging/configuration
set(CPACK_PACKAGE_VERSION ${SIMCOON_VERSION})

message(STATUS "Configuring ${PROJECT_NAME} ${SIMCOON_VERSION}")

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

# Options
option(SIMCOON_BUILD_TESTS "Build tests" ON)
option(SIMCOON_BUILD_PYTHON "Build Python bindings (for conda)" OFF)

# Detect conda build environment
if(DEFINED ENV{CONDA_BUILD})
  set(CONDA_BUILD TRUE)
  set(SIMCOON_BUILD_PYTHON ON)
  message(STATUS "Conda build detected - enabling Python bindings")
endif()

# Determine if we're building Python bindings (either via SKBUILD or SIMCOON_BUILD_PYTHON)
if(DEFINED SKBUILD OR SIMCOON_BUILD_PYTHON)
  set(BUILD_PYTHON_BINDINGS TRUE)
else()
  set(BUILD_PYTHON_BINDINGS FALSE)
endif()

# Helper function to set output directories for multi-config generators (MSVC)
function(set_output_directory target output_dir)
  set_target_properties(${target} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${output_dir}
    LIBRARY_OUTPUT_DIRECTORY ${output_dir}
    ARCHIVE_OUTPUT_DIRECTORY ${output_dir}
  )
  # For multi-config generators
  foreach(CONFIG_TYPE ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG_TYPE} CONFIG)
    set_target_properties(${target} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY_${CONFIG} ${output_dir}
      LIBRARY_OUTPUT_DIRECTORY_${CONFIG} ${output_dir}
      ARCHIVE_OUTPUT_DIRECTORY_${CONFIG} ${output_dir}
    )
  endforeach()
endfunction()

# Helper function to add test executable with proper configuration
# test_src should be relative to CMAKE_SOURCE_DIR (e.g., "test/Identification/Tidentification.cpp")
function(add_simcoon_test test_src)
  get_filename_component(testName ${test_src} NAME_WE)
  get_filename_component(testFolder ${test_src} DIRECTORY)

  # Strip base directory and intermediate grouping directories to match testBin structure:
  # test/Identification → Identification
  # test/Libraries/Continuum_mechanics → Continuum_mechanics
  # test/Umats/ELISO → Umats/ELISO
  # test_extern/Umats/UMABA → Umats/UMABA
  string(REGEX REPLACE "^test/" "" testFolder "${testFolder}")
  string(REGEX REPLACE "^test_extern/" "" testFolder "${testFolder}")
  string(REGEX REPLACE "^Libraries/" "" testFolder "${testFolder}")

  add_executable(${testName} ${CMAKE_SOURCE_DIR}/${test_src})
  target_link_libraries(${testName} PRIVATE simcoon GTest::gtest_main)

  set_output_directory(${testName} ${CMAKE_SOURCE_DIR}/testBin/${testFolder})

  if(MSVC)
    add_dependencies(${testName} simcoon)
    add_custom_command(TARGET ${testName} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_BINARY_DIR}/bin
      ${CMAKE_SOURCE_DIR}/testBin/${testFolder}
    )
  endif()

  # Create test output directory at build time (not configure time)
  add_custom_command(TARGET ${testName} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/testBin/${testFolder}
  )

  # Use generator expression to get correct executable path for multi-config generators
  add_test(NAME ${testFolder}/${testName}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testBin/${testFolder}
    COMMAND $<TARGET_FILE:${testName}>
  )
endfunction()

# Build type
if(NOT CMAKE_BUILD_TYPE)  # Debug by default
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel"
        FORCE)
endif()

message(STATUS "Build type ${CMAKE_BUILD_TYPE}")

# put our local cmake find scripts at the beginning of the cmake
# module search path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

#Setup CMake to run tests
enable_testing()

# Find dependencies

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# WORKAROUND for vcpkg BLAS/LAPACK debug/release library confusion on Windows
# vcpkg's FindBLAS/FindLAPACK modules (MODULE mode) don't respect CMAKE_BUILD_TYPE
# and may return debug libraries (from vcpkg/debug/lib/) even for Release builds.
# This causes ABI mismatches, crashes, and performance issues.
# Only apply fix for non-Debug builds to preserve proper debug library linking.
# TODO: Consider migrating to CONFIG mode packages (find_package(OpenBLAS CONFIG))
#       when conda/vcpkg compatibility issues are resolved.
if(WIN32 AND MSVC AND BLAS_FOUND AND LAPACK_FOUND)
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    string(REPLACE "/debug/lib/" "/lib/" BLAS_LIBRARIES_FIXED "${BLAS_LIBRARIES}")
    string(REPLACE "/debug/lib/" "/lib/" LAPACK_LIBRARIES_FIXED "${LAPACK_LIBRARIES}")

    # Only apply if we actually found debug paths (safety check)
    if(NOT BLAS_LIBRARIES_FIXED STREQUAL BLAS_LIBRARIES)
      set(BLAS_LIBRARIES ${BLAS_LIBRARIES_FIXED})
      set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES_FIXED})
      message(STATUS "Fixed BLAS/LAPACK: Using release libraries for ${CMAKE_BUILD_TYPE} build")
      message(STATUS "  BLAS: ${BLAS_LIBRARIES}")
      message(STATUS "  LAPACK: ${LAPACK_LIBRARIES}")
    endif()
  else()
    message(STATUS "Debug build: Using debug BLAS/LAPACK libraries as found")
  endif()
endif()

# Armadillo. Inside a conda environment, use the environment's own Armadillo (hence its
# BLAS/LAPACK and its single OpenMP runtime, conda-forge's llvm-openmp). Linking a
# system/Homebrew Armadillo into a conda Python process loads a second OpenMP runtime
# next to the one numpy/openblas/pytorch already carry: "OMP: Error #15" at import or
# random crashes in threaded kernels (KMP_DUPLICATE_LIB_OK only hides the guard).
# Install it with `conda install -c conda-forge armadillo`; -DArmadillo_ROOT=... still wins.
# Conda prefix: conda-build exposes the host prefix as PREFIX, an activated environment as
# CONDA_PREFIX (same precedence as the carma lookup below). `conda-meta/` is what makes a
# directory a conda prefix — PREFIX alone is a common Makefile variable (/usr/local ...).
set(SIMCOON_CONDA_PREFIX "")
if(DEFINED ENV{PREFIX} AND EXISTS "$ENV{PREFIX}/conda-meta")
  set(SIMCOON_CONDA_PREFIX "$ENV{PREFIX}")
elseif(DEFINED ENV{CONDA_PREFIX} AND EXISTS "$ENV{CONDA_PREFIX}/conda-meta")
  set(SIMCOON_CONDA_PREFIX "$ENV{CONDA_PREFIX}")
endif()
set(_simcoon_armadillo_auto OFF)
if(SIMCOON_CONDA_PREFIX AND NOT DEFINED Armadillo_ROOT AND NOT DEFINED ENV{Armadillo_ROOT})
  if(EXISTS "${SIMCOON_CONDA_PREFIX}/include/armadillo")
    set(Armadillo_ROOT "${SIMCOON_CONDA_PREFIX}")
    set(_simcoon_armadillo_auto ON)
    message(STATUS "Conda environment detected: Armadillo_ROOT = ${Armadillo_ROOT}")
  else()
    message(WARNING "Conda environment detected (${SIMCOON_CONDA_PREFIX}) but it has no Armadillo: "
                    "run `conda install -c conda-forge armadillo` so that simcoon shares the "
                    "environment's BLAS and OpenMP runtime (a non-conda Armadillo brings a second "
                    "OpenMP runtime into the Python process).")
  endif()
endif()
find_package(Armadillo 12.6 REQUIRED)
message(STATUS "Armadillo library: ${ARMADILLO_LIBRARIES}")
message(STATUS "Armadillo include: ${ARMADILLO_INCLUDE_DIRS}")
# literal search, not MATCHES: a prefix may contain regex metacharacters (py3.11, env(new))
string(FIND "${ARMADILLO_LIBRARIES};${ARMADILLO_INCLUDE_DIRS}" "${SIMCOON_CONDA_PREFIX}" _simcoon_arma_in_prefix)
if(_simcoon_armadillo_auto AND _simcoon_arma_in_prefix EQUAL -1)
  # a cached ARMADILLO_INCLUDE_DIR/ARMADILLO_LIBRARY from an earlier configure wins over the hint
  message(WARNING "Armadillo is NOT the one of the active conda environment (${SIMCOON_CONDA_PREFIX}): "
                  "${ARMADILLO_LIBRARIES}. Its BLAS will load a second OpenMP runtime next to the "
                  "environment's llvm-openmp (OMP Error #15 / crashes). Delete the stale build "
                  "cache (build/ or build/cp*) and configure again, or pass -DArmadillo_ROOT "
                  "explicitly if this is intended.")
endif()

# GTest is only needed when building tests
if(SIMCOON_BUILD_TESTS)
  find_package(GTest REQUIRED)
endif()

include(FetchContent)
FetchContent_Declare(
    dylib
    GIT_REPOSITORY "https://github.com/martin-olivier/dylib"
    GIT_TAG        "v3.0.1"
)

FetchContent_Declare(
    Fastor
    GIT_REPOSITORY "https://github.com/3MAH/Fastor.git"
    GIT_TAG        "b3499f6bea8bc49f92b5dcfc250da90e86fc5ebb"
)

#set(BUILD_SHARED_LIBS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
FetchContent_MakeAvailable(dylib)
set_target_properties(dylib PROPERTIES UNITY_BUILD OFF)

# dylib v3.0.1 relies on transitive includes that the libc++ shipped with the macOS 26.6
# SDK no longer pulls in: free() in demangle.cpp, errno in dylib.cpp, ntohl() in symbols.cpp.
# The runner image moved from Xcode 26.5 to 26.6 and the conda macOS build broke on it.
# Upstream fixed only the first one (martin-olivier/dylib 0a8573a9, unreleased), so force
# the headers here rather than moving the pin off a tagged release.
if(APPLE)
  # SHELL: is required: CMake de-duplicates repeated identical options, which would
  # collapse the three flags into "-include cstdlib cerrno arpa/inet.h" and make clang
  # treat the last two as input files ("no such file or directory: 'cerrno'").
  target_compile_options(dylib PRIVATE
    "SHELL:-include cstdlib" "SHELL:-include cerrno" "SHELL:-include arpa/inet.h")
endif()

# Fastor is header-only — just populate, no build step needed
FetchContent_Populate(Fastor)
message(STATUS "Fastor include dir: ${fastor_SOURCE_DIR}")

# Python dependencies (only if building Python bindings)
if(BUILD_PYTHON_BINDINGS)
  find_package(Python3 COMPONENTS Interpreter Development.Module NumPy REQUIRED)
  find_package(pybind11 REQUIRED)

  # pip/uv builds set SIMCOON_FETCH_CARMA (pyproject.toml) so every wheel and editable
  # install uses the pinned release fetched below, never a carma that happens to be on
  # the machine (e.g. another conda environment's). The conda recipe keeps its own.
  option(SIMCOON_FETCH_CARMA "Always fetch the pinned carma release instead of searching for one" OFF)

  # Try to find carma from conda/system first
  if(NOT SIMCOON_FETCH_CARMA)
    find_package(carma CONFIG QUIET)
  endif()

  if(carma_FOUND)
    message(STATUS "Found carma package (CONFIG mode)")
    # Get include directory from target for force-include
    get_target_property(carma_INCLUDE_DIR carma::carma INTERFACE_INCLUDE_DIRECTORIES)
    # Handle generator expressions - extract first real path
    if(carma_INCLUDE_DIR MATCHES "\\$<")
      set(carma_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
    endif()
    # Clear carma's INTERFACE_PRECOMPILE_HEADERS - conda packages have absolute paths
    # from the build machine that don't exist on the user's machine
    set_target_properties(carma::carma PROPERTIES INTERFACE_PRECOMPILE_HEADERS "")
  else()
    # If not found via CONFIG, try to find headers directly in conda prefix
    if(DEFINED ENV{PREFIX})
      set(CONDA_PREFIX $ENV{PREFIX})
    elseif(DEFINED ENV{CONDA_PREFIX})
      set(CONDA_PREFIX $ENV{CONDA_PREFIX})
    endif()

    if(NOT SIMCOON_FETCH_CARMA AND CONDA_PREFIX AND EXISTS "${CONDA_PREFIX}/include/carma")
      message(STATUS "Found carma headers in conda environment")
      set(carma_INCLUDE_DIR "${CONDA_PREFIX}/include")

      # Create carma::carma target from conda installation
      if(NOT TARGET carma::carma)
        add_library(carma::carma INTERFACE IMPORTED)
        set_target_properties(carma::carma PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${carma_INCLUDE_DIR}"
          INTERFACE_LINK_LIBRARIES Python3::NumPy
        )
      endif()
      set(carma_FOUND TRUE)
    else()
      # Fall back to fetching from GitHub
      message(STATUS "carma: fetching the pinned release v0.8.0 from GitHub")
      FetchContent_Declare(
        carma
        GIT_REPOSITORY https://github.com/RUrlus/carma.git
        GIT_TAG        v0.8.0
        # Headers only (no CMakeLists.txt there): carma's own project fetches a second
        # Armadillo (12.8.x) and puts it ahead of the one simcoon is built with, so _core
        # and libsimcoon would exchange arma objects across two Armadillo versions.
        SOURCE_SUBDIR  include
      )
      FetchContent_MakeAvailable(carma)

      # Get the source directory and set include path
      FetchContent_GetProperties(carma SOURCE_DIR carma_SOURCE_DIR)
      set(carma_INCLUDE_DIR "${carma_SOURCE_DIR}/include")

      # Ensure carma::carma target exists with proper include directories
      if(NOT TARGET carma::carma)
        add_library(carma::carma INTERFACE IMPORTED)
        set_target_properties(carma::carma PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${carma_INCLUDE_DIR}"
          INTERFACE_LINK_LIBRARIES Python3::NumPy
        )
      endif()
    endif()
  endif()

  message(STATUS "Python bindings enabled")
  message(STATUS "Python3_EXECUTABLE = ${Python3_EXECUTABLE}")
  message(STATUS "PYTHON_INCLUDE_DIRS = ${Python3_INCLUDE_DIRS}")
endif()

# Armadillo include directories will be added to targets as needed

# Set C++ standard globally
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Platform-specific configuration
if(MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

  # Set default output directories for multi-config generators
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  foreach(CONFIG_TYPE ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG_TYPE} CONFIG)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/bin")
  endforeach()
endif()

message(STATUS "CMAKE_SYSTEM_NAME          = ${CMAKE_SYSTEM_NAME}"         )
message(STATUS "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}"     )
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_CXX_STANDARD         = ${CMAKE_CXX_STANDARD}")
message(STATUS "CMAKE_COMPILER_IS_GNUCXX   = ${CMAKE_COMPILER_IS_GNUCXX}"  )
message(STATUS "BUILD_SHARED_LIBS          = ${BUILD_SHARED_LIBS}"         )
message(STATUS "DETECT_HDF5                = ${DETECT_HDF5}"               )
message(STATUS "ARMADILLO_INCLUDE_DIRS = ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "ARMADILLO_LIBRARIES = ${ARMADILLO_LIBRARIES}")

# Create simcoon shared library
add_library(simcoon SHARED)

# Set library version for ABI compatibility tracking
set_target_properties(simcoon PROPERTIES
  VERSION ${SIMCOON_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Add source files from subdirectory
add_subdirectory(src)

# Create alias for namespaced usage (simcoon::simcoon)
add_library(simcoon::simcoon ALIAS simcoon)

# Add include directories to simcoon target
target_include_directories(simcoon
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Add Armadillo as system headers to suppress warnings
target_include_directories(simcoon SYSTEM PUBLIC ${ARMADILLO_INCLUDE_DIRS})

# Add Fastor as system headers (header-only tensor library, build-time only)
target_include_directories(simcoon SYSTEM PUBLIC $<BUILD_INTERFACE:${fastor_SOURCE_DIR}>)


# Define restrict macro (MSVC and GCC/Clang handle restrict differently)
target_compile_definitions(simcoon PUBLIC restrict=)

# Platform-specific compile options
if(APPLE)
  target_compile_options(simcoon PUBLIC -stdlib=libc++)
  target_link_options(simcoon PUBLIC -stdlib=libc++)
endif()

# Optimization flags based on build type
if(NOT MSVC)
  target_compile_options(simcoon PUBLIC
    $<$<CONFIG:Debug>:-O0>
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:RelWithDebInfo>:-O2>
    $<$<CONFIG:MinSizeRel>:-Os>
  )
endif()

target_link_libraries(simcoon
    PUBLIC ${ARMADILLO_LIBRARIES}
    PRIVATE dylib
)

# OpenMP for batch tensor operations — NOT on Windows: libsimcoon is loaded into Python
# processes whose numpy/scipy stack carries another OpenMP runtime, and vcomp + iomp in one
# process deadlocks at the first parallel region (CI 6h hang in test_batch_apply_stress).
# Mirrors the python-builder threading policy: GCD on macOS, OpenMP on Linux, serial on Windows.
if(NOT WIN32)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(simcoon PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()

# On Windows with Python bindings, link libsimcoon to carma to use NumPy's memory allocator.
# This prevents heap corruption when Armadillo objects cross the libsimcoon/_core boundary.
# Only needed on Windows due to stricter heap validation.
# See: https://github.com/RUrlus/carma/issues/89
if(WIN32 AND BUILD_PYTHON_BINDINGS)
  target_link_libraries(simcoon PRIVATE carma::carma)

  # Force-include carma's cnalloc.h to ensure ARMA_ALIEN_MEM macros are defined
  # BEFORE any Armadillo headers.
  set(CARMA_CNALLOC_H "${carma_INCLUDE_DIR}/carma_bits/cnalloc.h")
  if(EXISTS "${CARMA_CNALLOC_H}")
    message(STATUS "Force-including carma cnalloc.h: ${CARMA_CNALLOC_H}")
    target_compile_options(simcoon PRIVATE "/FI${CARMA_CNALLOC_H}")
  else()
    message(WARNING "Could not find carma cnalloc.h at ${CARMA_CNALLOC_H}")
  endif()
endif()

##Testing
if(SIMCOON_BUILD_TESTS)
  # Test plugin libraries
  if(NOT MSVC)
    # Fortran is only needed for Abaqus UMAT plugin tests
    enable_language(Fortran)
    set_source_files_properties(testBin/Umats/UMABA/external/UMAT_ABAQUS_ELASTIC.f PROPERTIES Fortran_PREPROCESS OFF)
    add_library(umat_plugin_aba SHARED testBin/Umats/UMABA/external/UMAT_ABAQUS_ELASTIC.f testBin/Umats/UMABA/external/umat_plugin_aba.cpp)
    # Use default PREFIX (lib) and SUFFIX (.so on Linux, .dylib on macOS)
    # to match dylib::decorations::os_default() expectations
    set_target_properties(umat_plugin_aba PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMABA/external
    )
    target_link_libraries(umat_plugin_aba PRIVATE simcoon ${ARMADILLO_LIBRARIES})
  endif()

  add_library(umat_plugin_ext SHARED
    testBin/Umats/UMEXT/external/umat_plugin_ext.cpp
  )
  target_link_libraries(umat_plugin_ext PRIVATE simcoon)
  if(MSVC)
    set_target_properties(umat_plugin_ext PROPERTIES PREFIX "" SUFFIX ".dll")
    set_output_directory(umat_plugin_ext ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMEXT/external)
    add_dependencies(umat_plugin_ext simcoon)
    add_custom_command(TARGET umat_plugin_ext POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_BINARY_DIR}/bin
      ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMEXT/external
    )
  else()
    # Use default PREFIX (lib) and SUFFIX (.so on Linux, .dylib on macOS)
    # to match dylib::decorations::os_default() expectations
    set_target_properties(umat_plugin_ext PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMEXT/external
    )
  endif()

  # Add test executables from subdirectory
  add_subdirectory(test)
endif()

################################################################################
# INSTALL CONFIGURATION

message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMAKE_INSTALL_LIBDIR = ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "CMAKE_INSTALL_INCLUDEDIR = ${CMAKE_INSTALL_INCLUDEDIR}")
message(STATUS "CMAKE_INSTALL_BINDIR = ${CMAKE_INSTALL_BINDIR}")

# Find and collect runtime dependency DLLs on Windows (before Python bindings)
if(WIN32)
  # Find and install dependency DLLs (BLAS/LAPACK/runtime)
  set(DEPENDENCY_DLLS
      "openblas.dll"
      "liblapack.dll"
      "libgcc_s_seh-1.dll"
      "libgfortran-5.dll"
      "libquadmath-0.dll"
      "libwinpthread-1.dll"
  )

  # Search paths: library directories and compiler directory
  set(DLL_SEARCH_PATHS "")
  foreach(lib_path ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    get_filename_component(lib_dir ${lib_path} DIRECTORY)
    string(REPLACE "/lib" "/bin" bin_dir "${lib_dir}")
    list(APPEND DLL_SEARCH_PATHS ${bin_dir})
  endforeach()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_filename_component(COMPILER_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
    list(APPEND DLL_SEARCH_PATHS ${COMPILER_DIR})
  endif()

  list(REMOVE_DUPLICATES DLL_SEARCH_PATHS)

  # Collect dependency DLL paths for Python package
  if(BUILD_PYTHON_BINDINGS)
    set(PYTHON_RUNTIME_DLLS "" CACHE INTERNAL "Runtime DLLs needed by Python package")
  endif()

  # Try to find and store each DLL path
  foreach(dll_name ${DEPENDENCY_DLLS})
    set(dll_found FALSE)
    foreach(search_path ${DLL_SEARCH_PATHS})
      set(dll_path "${search_path}/${dll_name}")
      if(EXISTS "${dll_path}")
        # Store path for Python package (handled by python-builder CMakeLists)
        if(BUILD_PYTHON_BINDINGS)
          list(APPEND PYTHON_RUNTIME_DLLS "${dll_path}")
          set(PYTHON_RUNTIME_DLLS ${PYTHON_RUNTIME_DLLS} CACHE INTERNAL "Runtime DLLs needed by Python package")
        endif()
        message(STATUS "Found dependency DLL: ${dll_path}")
        set(dll_found TRUE)
        break()
      endif()
    endforeach()

    if(NOT dll_found)
      message(WARNING "Dependency DLL not found: ${dll_name}. "
              "Python package may fail at runtime. "
              "Search paths: ${DLL_SEARCH_PATHS}")
    endif()
  endforeach()
endif()

################################################################################
# PYTHON BINDINGS
if(BUILD_PYTHON_BINDINGS)
  add_subdirectory(simcoon-python-builder)
endif()

################################################################################
# INSTALL TARGETS
# Skip C++ library installation when building Python wheels (SKBUILD)
# The Python builder handles installing libsimcoon.so to the wheel package

if(NOT DEFINED SKBUILD)

# Install simcoon library with export
install(TARGETS simcoon
  EXPORT simcoonTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT Runtime
    NAMELINK_COMPONENT Development
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT Development
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT Runtime
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install headers
install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development
)

# Install the export set
install(EXPORT simcoonTargets
  FILE simcoonTargets.cmake
  NAMESPACE simcoon::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/simcoon
  COMPONENT Development
)

# Generate and install package config files
include(CMakePackageConfigHelpers)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/simcoonConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/simcoonConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/simcoon
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/simcoonConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/simcoonConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/simcoonConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/simcoon
  COMPONENT Development
)

# Install runtime dependencies on Windows
if(WIN32)
  # Install all DLLs from build/bin
  install(DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
          DESTINATION ${CMAKE_INSTALL_BINDIR}
          COMPONENT Runtime
          FILES_MATCHING PATTERN "*.dll")

  # Install dependency DLLs
  if(PYTHON_RUNTIME_DLLS)
    install(FILES ${PYTHON_RUNTIME_DLLS}
            DESTINATION ${CMAKE_INSTALL_BINDIR}
            COMPONENT Runtime)
  endif()
endif()

endif() # NOT DEFINED SKBUILD
