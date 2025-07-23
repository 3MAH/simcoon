# Main CMakeLists.txt

cmake_minimum_required(VERSION 3.21)
project(simcoon LANGUAGES C CXX)

# Project version from version file
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/python-setup/simcoon/__version__.py" VERSION_FILE)
string(REGEX MATCH "__version__ = \"([^\"]+)\"" _ "${VERSION_FILE}")
set(PROJECT_VERSION "${CMAKE_MATCH_1}")
message(STATUS "Configuring ${PROJECT_NAME} ${PROJECT_VERSION}")

# Options
option(DEBUG "Compile with debugging information" OFF)
option(BUILD_TESTS "Build the tests" ON)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# For multi-config generators (Visual Studio)
foreach(CONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG} CONFIG_UPPER)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG_UPPER} ${CMAKE_BINARY_DIR}/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG_UPPER} ${CMAKE_BINARY_DIR}/lib)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG_UPPER} ${CMAKE_BINARY_DIR}/lib)
endforeach()

# Platform-specific flags
if(MSVC)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict=")
else()
    enable_language(Fortran)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -Drestrict=")
    if(NOT DEBUG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    endif()
    if(APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -fpermissive")
    endif()
endif()

# Dependencies
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Armadillo 12.6 REQUIRED)

# Boost detection - Linux-specific fix for manylinux environments
if(UNIX AND NOT APPLE)
    # Linux: Use manual detection for manylinux compatibility
    message(STATUS "Linux detected: Using manual Boost detection for manylinux compatibility")
    
    # Check if we're in a manylinux environment (has the specific library paths)
    if(EXISTS "/usr/lib64/libboost_filesystem.so" AND EXISTS "/usr/lib64/libboost_system.so")
        # Direct linking for manylinux
        set(Boost_FOUND TRUE)
        set(Boost_INCLUDE_DIRS /usr/include)
        set(Boost_LIBRARIES /usr/lib64/libboost_filesystem.so /usr/lib64/libboost_system.so)
        message(STATUS "Using direct Boost linking for manylinux:")
        message(STATUS "  Include dirs: ${Boost_INCLUDE_DIRS}")
        message(STATUS "  Libraries: ${Boost_LIBRARIES}")
    else()
        # Standard Linux detection
        find_package(Boost 1.57.0 COMPONENTS filesystem system)
        if(NOT Boost_FOUND)
            # Fallback manual detection
            find_path(Boost_INCLUDE_DIR 
                NAMES boost/version.hpp
                PATHS /usr/include /usr/local/include)
            find_library(Boost_FILESYSTEM_LIBRARY 
                NAMES boost_filesystem
                PATHS /usr/lib64 /usr/lib /usr/local/lib)
            find_library(Boost_SYSTEM_LIBRARY 
                NAMES boost_system  
                PATHS /usr/lib64 /usr/lib /usr/local/lib)
            
            if(Boost_INCLUDE_DIR AND Boost_FILESYSTEM_LIBRARY AND Boost_SYSTEM_LIBRARY)
                set(Boost_FOUND TRUE)
                set(Boost_LIBRARIES ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
                set(Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIR})
                message(STATUS "Manual Boost detection successful")
            else()
                message(FATAL_ERROR "Could not find Boost libraries on Linux")
            endif()
        endif()
    endif()
else()
    # Windows and macOS: Use standard detection (already working)
    find_package(Boost 1.57.0)
endif()

# Include directories
include_directories(include)
include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
if(Boost_FOUND)
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
endif()

# Main library
file(GLOB_RECURSE source_files src*/* include/*)
add_library(simcoon SHARED ${source_files})

# Link libraries - include Boost on Linux, optional on other platforms
if(UNIX AND NOT APPLE)
    # Linux: Always link Boost libraries
    target_link_libraries(simcoon ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES})
else()
    # Windows/macOS: Original logic
    target_link_libraries(simcoon ${ARMADILLO_LIBRARIES})
    if(Boost_FOUND)
        target_link_libraries(simcoon ${Boost_LIBRARIES})
    endif()
endif()

if(MSVC)
    target_link_libraries(simcoon ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    if(BUILD_PYTHON_BINDINGS)
        target_link_libraries(simcoon carma::carma)
    endif()
endif()

# Executables with proper RPATH configuration
set(EXECUTABLES solver identification L_eff Elastic_props ODF PDF)
foreach(exe ${EXECUTABLES})
    add_executable(${exe} software/${exe}.cpp)
    
    # Link libraries - include Boost on Linux
    if(UNIX AND NOT APPLE)
        # Linux: Always link Boost libraries
        target_link_libraries(${exe} simcoon ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES})
    else()
        # Windows/macOS: Original logic
        target_link_libraries(${exe} simcoon ${ARMADILLO_LIBRARIES})
        if(Boost_FOUND)
            target_link_libraries(${exe} ${Boost_LIBRARIES})
        endif()
    endif()
    
    if(NOT WIN32)
        target_link_libraries(${exe} -ldl)
    endif()
    
    # Configure RPATH for all platforms
    if(APPLE)
        set_target_properties(${exe} PROPERTIES 
            INSTALL_RPATH "@loader_path"
            BUILD_WITH_INSTALL_RPATH TRUE
            # Set deployment target for compatibility
            OSX_DEPLOYMENT_TARGET "${CMAKE_OSX_DEPLOYMENT_TARGET}")
    elseif(UNIX)
        set_target_properties(${exe} PROPERTIES 
            INSTALL_RPATH "$ORIGIN"
            BUILD_WITH_INSTALL_RPATH TRUE)
    endif()
endforeach()

# Testing
if(BUILD_TESTS)
    enable_testing()
    
    # GTest
    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.17.0
        )
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
        add_library(GTest::gtest ALIAS gtest)
        add_library(GTest::gtest_main ALIAS gtest_main)
    endif()
    
    # Test plugins
    if(NOT MSVC)
        add_library(umat_plugin_aba SHARED 
            testBin/Umats/UMABA/external/UMAT_ABAQUS_ELASTIC.for 
            testBin/Umats/UMABA/external/umat_plugin_aba.cpp)
        set_target_properties(umat_plugin_aba PROPERTIES 
            PREFIX "" SUFFIX "" 
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMABA/external)
        
        # Link with Boost on Linux
        if(UNIX AND NOT APPLE)
            target_link_libraries(umat_plugin_aba simcoon ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES})
        else()
            target_link_libraries(umat_plugin_aba simcoon ${ARMADILLO_LIBRARIES})
        endif()
    endif()
    
    add_library(umat_plugin_ext SHARED testBin/Umats/UMEXT/external/umat_plugin_ext.cpp)
    
    # Link with Boost on Linux
    if(UNIX AND NOT APPLE)
        target_link_libraries(umat_plugin_ext simcoon ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES})
    else()
        target_link_libraries(umat_plugin_ext simcoon ${ARMADILLO_LIBRARIES})
    endif()
    
    if(MSVC)
        set_target_properties(umat_plugin_ext PROPERTIES 
            PREFIX "" 
            SUFFIX ".dll"
            RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMEXT/external
            RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMEXT/external)
    else()
        set_target_properties(umat_plugin_ext PROPERTIES 
            PREFIX "" SUFFIX "" 
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMEXT/external)
    endif()
    
    # Test executables
    file(GLOB_RECURSE TEST_SRCS test/*.cpp)
    foreach(testSrc ${TEST_SRCS})
        get_filename_component(testName ${testSrc} NAME_WE)
        file(RELATIVE_PATH testRel ${CMAKE_CURRENT_SOURCE_DIR}/test ${testSrc})
        get_filename_component(testFolder ${testRel} DIRECTORY)
        
        add_executable(${testName} ${testSrc})
        
        # Link with Boost on Linux
        if(UNIX AND NOT APPLE)
            target_link_libraries(${testName} simcoon ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} GTest::gtest GTest::gtest_main)
        else()
            target_link_libraries(${testName} simcoon ${ARMADILLO_LIBRARIES} GTest::gtest GTest::gtest_main)
        endif()
        
        set_target_properties(${testName} PROPERTIES 
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
            RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
            RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder})
        
        if(MSVC)
            add_custom_command(TARGET ${testName} POST_BUILD 
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE_DIR:simcoon>/simcoon.dll  
                ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}/)
        endif()
        
        add_test(NAME ${testFolder}/${testName} 
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}/${testName})
    endforeach()
endif()

# Python bindings
if(BUILD_PYTHON_BINDINGS)
    message(STATUS "Building Python bindings...")
    
    # Find Python with all required components
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
    
    # Find pybind11 first
    find_package(pybind11 CONFIG QUIET)
    if(NOT pybind11_FOUND)
        # Try to find pybind11 through Python
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
            OUTPUT_VARIABLE pybind11_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE pybind11_RESULT
        )
        if(pybind11_RESULT EQUAL 0)
            # Convert Windows backslashes to forward slashes
            string(REPLACE "\\" "/" pybind11_DIR "${pybind11_DIR}")
            find_package(pybind11 CONFIG REQUIRED PATHS "${pybind11_DIR}")
        else()
            message(FATAL_ERROR "Could not find pybind11. Please install it: pip install pybind11")
        endif()
    endif()
    
    # Now find/fetch carma - it will use the pybind11 we already found
    set(CARMA_USE_EXTERNAL_PYBIND11 ON CACHE BOOL "Use external pybind11" FORCE)
    find_package(carma CONFIG QUIET)
    if(NOT carma_FOUND)
        include(FetchContent)
        FetchContent_Declare(carma
            GIT_REPOSITORY https://github.com/RUrlus/carma.git
            GIT_TAG v0.8.0)
        FetchContent_MakeAvailable(carma)
    endif()
    
    # Add Python builder
    add_subdirectory(simcoon-python-builder)
endif()

# Installation
include(GNUInstallDirs)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS simcoon)

if(BUILD_PYTHON_BINDINGS)
    if(WIN32)
        install(TARGETS simmit simcoon ${EXECUTABLES}
            DESTINATION simcoon 
            COMPONENT python)
        
        # Only search in bin directory to avoid conflicts
        set(REQUIRED_DLLS
            openblas.dll liblapack.dll libgcc_s_seh-1.dll 
            libgfortran-5.dll libquadmath-0.dll libwinpthread-1.dll)
        
        foreach(dll ${REQUIRED_DLLS})
            install(FILES ${CMAKE_BINARY_DIR}/bin/${dll}
                DESTINATION simcoon
                COMPONENT python
                OPTIONAL)
        endforeach()
    else()
        install(TARGETS simmit simcoon ${EXECUTABLES} 
            DESTINATION simcoon 
            COMPONENT python)
    endif()
    
    install(DIRECTORY python-setup/simcoon/ 
        DESTINATION simcoon 
        COMPONENT python
        FILES_MATCHING PATTERN "*.py")
endif()