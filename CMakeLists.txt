cmake_minimum_required(VERSION 3.21)
project(simcoon LANGUAGES C CXX)

file(READ "python-setup/simcoon/__version__.py" VERSION_FILE)
string(REGEX MATCH "__version__ = \"([^\"]+)\"" _ "${VERSION_FILE}")
set(PROJECT_VERSION "${CMAKE_MATCH_1}")
message(STATUS "Configuring ${PROJECT_NAME} ${PROJECT_VERSION}")

option(BUILD_TESTS "Build the tests" ON)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Armadillo 12.6 REQUIRED)
find_package(Boost 1.57.0)

file(GLOB_RECURSE SOURCE_FILES 
    "src/*/*.cpp" 
    "src/*/*.c"
    "include/*.hpp"
    "include/*.h"
)

add_library(simcoon SHARED ${SOURCE_FILES})

target_include_directories(simcoon 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${ARMADILLO_INCLUDE_DIRS}
)

target_link_libraries(simcoon 
    PUBLIC ${ARMADILLO_LIBRARIES}
    PRIVATE 
        ${BLAS_LIBRARIES} 
        ${LAPACK_LIBRARIES}
)

if(MSVC)
    target_compile_definitions(simcoon PRIVATE restrict=)
    set_target_properties(simcoon PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
    enable_language(Fortran)
    target_compile_definitions(simcoon PRIVATE restrict=)
    target_compile_options(simcoon PRIVATE 
        $<$<CONFIG:Release>:-O3>
        $<$<PLATFORM_ID:Darwin>:-stdlib=libc++ -fpermissive>
    )
endif()

if(Boost_FOUND)
    target_include_directories(simcoon PRIVATE ${Boost_INCLUDE_DIRS})
endif()

function(add_simcoon_executable name)
    add_executable(${name} software/${name}.cpp)
    target_link_libraries(${name} PRIVATE simcoon)
    if(UNIX)
        target_link_libraries(${name} PRIVATE dl)
    endif()
endfunction()

set(EXECUTABLES solver identification L_eff Elastic_props ODF PDF)
foreach(exe ${EXECUTABLES})
    add_simcoon_executable(${exe})
endforeach()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(simcoon-python-builder)
endif()

include(GNUInstallDirs)
install(TARGETS simcoon ${EXECUTABLES}
    EXPORT simcoonTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/ 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT simcoonTargets
    FILE simcoonTargets.cmake
    NAMESPACE simcoon::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/simcoon
)