cmake_minimum_required(VERSION 3.21)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
        OUTPUT_VARIABLE pybind11_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE pybind11_RESULT
    )
    if(pybind11_RESULT EQUAL 0)
        string(REPLACE "\\" "/" pybind11_DIR "${pybind11_DIR}")
        find_package(pybind11 CONFIG REQUIRED PATHS "${pybind11_DIR}")
    else()
        message(FATAL_ERROR "pybind11 not found. Install with: pip install pybind11")
    endif()
endif()

find_package(carma CONFIG QUIET)
if(NOT carma_FOUND)
    include(FetchContent)
    FetchContent_Declare(carma
        GIT_REPOSITORY https://github.com/RUrlus/carma.git
        GIT_TAG v0.8.0
    )
    set(CARMA_USE_EXTERNAL_PYBIND11 ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(carma)
endif()

file(GLOB_RECURSE SIMMIT_SOURCES "src/python_wrappers/*.cpp")
if(NOT SIMMIT_SOURCES)
    message(FATAL_ERROR "No Python wrapper sources found")
endif()

pybind11_add_module(simmit ${SIMMIT_SOURCES})

target_link_libraries(simmit PRIVATE 
    simcoon 
    carma::carma
)

target_include_directories(simmit PRIVATE 
    include
    ${CMAKE_SOURCE_DIR}/include
)

if(MSVC)
    add_custom_command(TARGET simmit POST_BUILD 
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:simcoon>
        $<TARGET_FILE_DIR:simmit>
    )
else()
    set_target_properties(simmit PROPERTIES 
        INSTALL_RPATH "$<IF:$<PLATFORM_ID:Darwin>,@loader_path,$ORIGIN>"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()

if(BUILD_TESTS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/Carma/test_carma.cpp")
    pybind11_add_module(test_carma test/Carma/test_carma.cpp)
    target_link_libraries(test_carma PRIVATE carma::carma)
    
    if(NOT MSVC)
        set_target_properties(test_carma PROPERTIES
            PREFIX ""
            SUFFIX ".so"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test/Carma"
        )
    endif()
    
    foreach(test_script run_test_roundtrip.py run_test_arr_to_mat.py)
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/Carma/${test_script}")
            string(REPLACE ".py" "" test_name ${test_script})
            add_test(NAME "carma_${test_name}"
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test/Carma"
                COMMAND ${Python3_EXECUTABLE} ${test_script}
            )
        endif()
    endforeach()
endif()

if(BUILD_PYTHON_BINDINGS)
    install(TARGETS simmit simcoon ${EXECUTABLES}
        DESTINATION simcoon 
        COMPONENT python
    )
    
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/python-setup/simcoon/ 
        DESTINATION simcoon 
        COMPONENT python
        FILES_MATCHING PATTERN "*.py"
    )
    
    if(WIN32)
        set(REQUIRED_DLLS
            openblas.dll liblapack.dll libgcc_s_seh-1.dll 
            libgfortran-5.dll libquadmath-0.dll libwinpthread-1.dll
        )
        
        foreach(dll ${REQUIRED_DLLS})
            install(FILES "${CMAKE_BINARY_DIR}/bin/${dll}"
                DESTINATION simcoon
                COMPONENT python
                OPTIONAL
            )
        endforeach()
    endif()
endif()