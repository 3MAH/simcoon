cmake_minimum_required(VERSION 3.21)

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
endif()

function(add_test_plugin name sources)
    add_library(${name} SHARED ${sources})
    target_link_libraries(${name} PRIVATE simcoon)
    
    set_target_properties(${name} PROPERTIES 
        PREFIX "" 
        SUFFIX "$<IF:$<PLATFORM_ID:Windows>,.dll,>"
    )
    
    if(name STREQUAL "umat_plugin_aba")
        set_target_properties(${name} PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../testBin/Umats/UMABA/external"
        )
    elseif(name STREQUAL "umat_plugin_ext")
        if(MSVC)
            set_target_properties(${name} PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../testBin/Umats/UMEXT/external"
            )
        else()
            set_target_properties(${name} PROPERTIES
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../testBin/Umats/UMEXT/external"
            )
        endif()
    endif()
endfunction()

if(NOT MSVC)
    add_test_plugin(umat_plugin_aba 
        "../testBin/Umats/UMABA/external/UMAT_ABAQUS_ELASTIC.for;../testBin/Umats/UMABA/external/umat_plugin_aba.cpp"
    )
endif()

add_test_plugin(umat_plugin_ext 
    "../testBin/Umats/UMEXT/external/umat_plugin_ext.cpp"
)

function(add_simcoon_test test_source)
    get_filename_component(test_name ${test_source} NAME_WE)
    file(RELATIVE_PATH test_rel "${CMAKE_CURRENT_SOURCE_DIR}/../test" ${test_source})
    get_filename_component(test_folder ${test_rel} DIRECTORY)
    
    add_executable(${test_name} ${test_source})
    target_link_libraries(${test_name} PRIVATE 
        simcoon 
        GTest::gtest 
        GTest::gtest_main
    )
    
    set(test_output_dir "${CMAKE_CURRENT_SOURCE_DIR}/../testBin/${test_folder}")
    set_target_properties(${test_name} PROPERTIES 
        RUNTIME_OUTPUT_DIRECTORY "${test_output_dir}"
    )
    
    if(MSVC)
        add_custom_command(TARGET ${test_name} POST_BUILD 
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:simcoon>  
            "${test_output_dir}/"
        )
    endif()
    
    add_test(NAME "${test_folder}/${test_name}" 
        WORKING_DIRECTORY "${test_output_dir}"
        COMMAND "${test_output_dir}/${test_name}"
    )
endfunction()

file(GLOB_RECURSE TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/../test/*.cpp")
foreach(test_src ${TEST_SOURCES})
    add_simcoon_test(${test_src})
endforeach()