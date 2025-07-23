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
find_package(Boost 1.57.0)

# Include directories
include_directories(include)
include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
if(Boost_FOUND)
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
endif()

# Main library
file(GLOB_RECURSE source_files src*/* include/*)
add_library(simcoon SHARED ${source_files})
target_link_libraries(simcoon ${ARMADILLO_LIBRARIES})
if(MSVC)
    target_link_libraries(simcoon ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    if(BUILD_PYTHON_BINDINGS)
        target_link_libraries(simcoon carma::carma)
    endif()
endif()

# Executables
set(EXECUTABLES solver identification L_eff Elastic_props ODF PDF)
foreach(exe ${EXECUTABLES})
    add_executable(${exe} software/${exe}.cpp)
    target_link_libraries(${exe} simcoon ${ARMADILLO_LIBRARIES})
    if(NOT WIN32)
        target_link_libraries(${exe} -ldl)
    endif()
    
    # Add RPATH configuration for macOS
    if(APPLE)
        set_target_properties(${exe} PROPERTIES 
            INSTALL_RPATH "@loader_path"
            BUILD_WITH_INSTALL_RPATH TRUE)
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
        target_link_libraries(umat_plugin_aba simcoon ${ARMADILLO_LIBRARIES})
    endif()
    
    add_library(umat_plugin_ext SHARED testBin/Umats/UMEXT/external/umat_plugin_ext.cpp)
    target_link_libraries(umat_plugin_ext simcoon ${ARMADILLO_LIBRARIES})
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
        target_link_libraries(${testName} simcoon ${ARMADILLO_LIBRARIES} GTest::gtest GTest::gtest_main)
        
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


# if(BUILD_PYTHON_BINDINGS)
#     if(WIN32)
#         # Install the main Python module and simcoon.dll
#         install(TARGETS simmit 
#             RUNTIME DESTINATION simcoon 
#             LIBRARY DESTINATION simcoon
#             COMPONENT python)
#         install(TARGETS simcoon
#             RUNTIME DESTINATION simcoon 
#             COMPONENT python)
#         install(TARGETS ${EXECUTABLES}
#             RUNTIME DESTINATION simcoon 
#             COMPONENT python)
            
#         # Install required runtime DLLs using install-time code execution
#         install(CODE "
#             set(REQUIRED_DLL_NAMES
#                 \"openblas.dll\"
#                 \"liblapack.dll\"
#                 \"libgcc_s_seh-1.dll\"
#                 \"libgfortran-5.dll\"
#                 \"libquadmath-0.dll\"
#                 \"libwinpthread-1.dll\"
#             )
            
#             # Search paths - these will be evaluated at install time
#             set(SEARCH_PATHS
#                 \"${CMAKE_BINARY_DIR}/bin\"
#                 \"${CMAKE_RUNTIME_OUTPUT_DIRECTORY}\"
#             )
            
#             message(STATUS \"Installing required DLLs at install time...\")
            
#             foreach(DLL_NAME \${REQUIRED_DLL_NAMES})
#                 set(DLL_FOUND FALSE)
#                 foreach(SEARCH_PATH \${SEARCH_PATHS})
#                     # Use file(GLOB) to find the DLL at install time
#                     file(GLOB DLL_MATCHES \"\${SEARCH_PATH}/\${DLL_NAME}\")
#                     if(DLL_MATCHES)
#                         list(GET DLL_MATCHES 0 DLL_PATH)
#                         message(STATUS \"Installing \${DLL_NAME} from \${DLL_PATH}\")
#                         file(INSTALL \"\${DLL_PATH}\" 
#                              DESTINATION \"\${CMAKE_INSTALL_PREFIX}/simcoon\"
#                              USE_SOURCE_PERMISSIONS)
#                         set(DLL_FOUND TRUE)
#                         break()
#                     endif()
#                 endforeach()
                
#                 if(NOT DLL_FOUND)
#                     message(WARNING \"Required DLL not found at install time: \${DLL_NAME}\")
#                     message(STATUS \"Searched in: \${SEARCH_PATHS}\")
#                     # Debug: List what files are actually in the bin directory
#                     file(GLOB BIN_CONTENTS \"${CMAKE_BINARY_DIR}/bin/*.dll\")
#                     message(STATUS \"DLLs found in bin directory:\")
#                     foreach(dll_file \${BIN_CONTENTS})
#                         message(STATUS \"  \${dll_file}\")
#                     endforeach()
#                 endif()
#             endforeach()
#         " COMPONENT python)
            
#     else()
#         install(TARGETS simmit simcoon ${EXECUTABLES} 
#             DESTINATION simcoon 
#             COMPONENT python)
#     endif()
    
#     install(DIRECTORY python-setup/simcoon/ 
#         DESTINATION simcoon 
#         COMPONENT python
#         FILES_MATCHING PATTERN "*.py")
# endif()