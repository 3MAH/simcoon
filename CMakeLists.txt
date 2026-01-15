# CMakeLists.txt

cmake_minimum_required(VERSION 3.15)

project(simcoon LANGUAGES C CXX)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)

# Read version from Python package
file(READ "python-setup/simcoon/__version__.py" VERSION_FILE)
string(REGEX MATCH "__version__ = \"([0-9]+\\.[0-9]+\\.[0-9]+)\"" _ ${VERSION_FILE})
set(PROJECT_VERSION ${CMAKE_MATCH_1})

# Extract major version for SOVERSION
string(REGEX MATCH "^([0-9]+)" _ ${PROJECT_VERSION})
set(PROJECT_VERSION_MAJOR ${CMAKE_MATCH_1})

# Print the project version
message(STATUS "Project version: ${PROJECT_VERSION}")

# Optionally use this version in packaging/configuration
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})

message(STATUS "Configuring ${PROJECT_NAME} ${PROJECT_VERSION}")

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

# Options
option(SIMCOON_BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(SIMCOON_BUILD_TESTS "Build tests" ON)

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

  # When Python bindings are enabled, libsimcoon links to carma which requires Python symbols.
  # Link tests to Python so they can resolve these symbols at runtime.
  if(SIMCOON_BUILD_PYTHON_BINDINGS)
    target_link_libraries(${testName} PRIVATE Python3::Python)
  endif()

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
if(MSVC)
  set(BLA_PKGCONFIG_BLAS True)
endif()

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

find_package(Armadillo 12.6 REQUIRED)

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

#set(BUILD_SHARED_LIBS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
FetchContent_MakeAvailable(dylib)
set_target_properties(dylib PROPERTIES UNITY_BUILD OFF)

# Python dependencies (only if building Python bindings)
if(SIMCOON_BUILD_PYTHON_BINDINGS)
  find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
  find_package(pybind11 REQUIRED)

  # Try to find carma from conda/system first
  find_package(carma CONFIG QUIET)

  if(carma_FOUND)
    message(STATUS "Found carma package (CONFIG mode)")
  else()
    # If not found via CONFIG, try to find headers directly in conda prefix
    if(DEFINED ENV{PREFIX})
      set(CONDA_PREFIX $ENV{PREFIX})
    elseif(DEFINED ENV{CONDA_PREFIX})
      set(CONDA_PREFIX $ENV{CONDA_PREFIX})
    endif()

    if(CONDA_PREFIX AND EXISTS "${CONDA_PREFIX}/include/carma")
      message(STATUS "Found carma headers in conda environment")
      set(carma_INCLUDE_DIR "${CONDA_PREFIX}/include")

      # Create carma::carma target from conda installation
      if(NOT TARGET carma::carma)
        add_library(carma::carma INTERFACE IMPORTED)
        set_target_properties(carma::carma PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${carma_INCLUDE_DIR}"
        )
      endif()
      set(carma_FOUND TRUE)
    else()
      # Fall back to fetching from GitHub
      message(STATUS "carma not found, fetching from GitHub")
      FetchContent_Declare(
        carma
        GIT_REPOSITORY https://github.com/RUrlus/carma.git
        GIT_TAG        v0.8.0
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
else()
  # Unix-like systems (Linux, macOS)
  enable_language(Fortran)
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
  VERSION ${PROJECT_VERSION}
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

# When building Python bindings, link libsimcoon to carma::carma
# This ensures libsimcoon uses NumPy's memory allocator (via ARMA_ALIEN_MEM_ALLOC_FUNCTION)
# which prevents heap corruption on Windows when Armadillo objects cross the
# libsimcoon/simmit boundary. See: https://github.com/RUrlus/carma/issues/89
if(SIMCOON_BUILD_PYTHON_BINDINGS)
  target_link_libraries(simcoon PUBLIC carma::carma)
endif()

# Compile public executables (with main functions)
set(EXECUTABLES solver identification L_eff Elastic_props ODF PDF)
foreach(exe_name ${EXECUTABLES})
  add_executable(${exe_name} software/${exe_name}.cpp)
  target_link_libraries(${exe_name} PRIVATE simcoon ${CMAKE_DL_LIBS})
endforeach()

##Testing
if(SIMCOON_BUILD_TESTS)
  # Test plugin libraries
  if(NOT MSVC)
    add_library(umat_plugin_aba SHARED testBin/Umats/UMABA/external/UMAT_ABAQUS_ELASTIC.for testBin/Umats/UMABA/external/umat_plugin_aba.cpp)

    set_target_properties(umat_plugin_aba PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMABA/external
    )

    target_link_libraries(umat_plugin_aba PRIVATE simcoon)
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
  if(SIMCOON_BUILD_PYTHON_BINDINGS)
    set(PYTHON_RUNTIME_DLLS "" CACHE INTERNAL "Runtime DLLs needed by Python package")
  endif()

  # Try to find and store each DLL path
  foreach(dll_name ${DEPENDENCY_DLLS})
    set(dll_found FALSE)
    foreach(search_path ${DLL_SEARCH_PATHS})
      set(dll_path "${search_path}/${dll_name}")
      if(EXISTS "${dll_path}")
        # Store path for Python package (handled by python-builder CMakeLists)
        if(SIMCOON_BUILD_PYTHON_BINDINGS)
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
if(SIMCOON_BUILD_PYTHON_BINDINGS)
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
