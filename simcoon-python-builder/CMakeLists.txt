# simcoon-python-builder/CMakeLists.txt

cmake_minimum_required(VERSION 3.15)
project(simcoon-python-builder)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Collect source files for the Python module
file(GLOB_RECURSE SIMMIT_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/python_wrappers/*.cpp
)

# Check if we have source files
if(NOT SIMMIT_SOURCES)
    message(FATAL_ERROR "No source files found in ${CMAKE_CURRENT_SOURCE_DIR}/src/python_wrappers/")
endif()

# Python module
pybind11_add_module(simmit ${SIMMIT_SOURCES})

# Platform-specific extension
if(MSVC)
    set_target_properties(simmit PROPERTIES SUFFIX ".pyd")
endif()

# Link dependencies
target_link_libraries(simmit PRIVATE
    simcoon
    carma::carma
    ${ARMADILLO_LIBRARIES}
)

# Include directories for the Python module
target_include_directories(simmit PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include
)

# RPATH settings
if(APPLE)
    set_target_properties(simmit PROPERTIES 
        INSTALL_RPATH "@loader_path"
        BUILD_WITH_INSTALL_RPATH TRUE)
elseif(UNIX)
    set_target_properties(simmit PROPERTIES 
        INSTALL_RPATH "$ORIGIN"
        BUILD_WITH_INSTALL_RPATH TRUE)
endif()

# Testing
if(BUILD_TESTS)
    # Check if test file exists
    set(TEST_CARMA_FILE ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma/test_carma.cpp)
    if(EXISTS ${TEST_CARMA_FILE})
        # Test module
        pybind11_add_module(test_carma ${TEST_CARMA_FILE})
        
        if(MSVC)
            set_target_properties(test_carma PROPERTIES SUFFIX ".pyd")
        else()
            set_target_properties(test_carma PROPERTIES
                PREFIX ""
                SUFFIX ".so"
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test/Carma"
            )
        endif()
        
        target_link_libraries(test_carma PRIVATE
            carma::carma
            ${ARMADILLO_LIBRARIES}
        )
        
        # Add tests only if the Python test scripts exist
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma/run_test_roundtrip.py)
            add_test(NAME test_carma_roundtrip
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma
                COMMAND ${Python3_EXECUTABLE} run_test_roundtrip.py
            )
        endif()
        
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma/run_test_arr_to_mat.py)
            add_test(NAME test_carma_arr_to_mat
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma
                COMMAND ${Python3_EXECUTABLE} run_test_arr_to_mat.py
            )
        endif()
    else()
        message(STATUS "Test file ${TEST_CARMA_FILE} not found. Skipping carma tests.")
    endif()
endif()