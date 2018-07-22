# CMakeLists.txt

cmake_minimum_required(VERSION 2.6)
project(simcoon)

# The version number
set(VERSION_MAJOR 1)
set(VERSION_MINOR 9)
set(VERSION_PATCH 2)

if(VERSION_PATCH MATCHES "0")
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}")
else()
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()

message(STATUS "Configuring ${PROJECT_NAME} ${VERSION_NUMBER}")

# We default to debugging mode for developers.
option(DEBUG "Compile with debugging information" OFF)

# configure a header file to pass some CMake settings
# to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/configuration.hpp.in"
#  "${PROJECT_BINARY_DIR}/configuration.hpp"
#)

# put our local cmake find scripts at the beginning of the cmake
# module search path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

#Setup CMake to run tests
enable_testing()

#Set executable files and library files
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

#Inclusion of Armadillo
include(FindArmadillo)
find_package(Armadillo 5.2 REQUIRED)
include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})

# Find FTensor
find_package(FTensor)
include_directories(SYSTEM ${FTENSOR_INCLUDE_DIRS})

# Find CGAL
find_package(CGAL REQUIRED COMPONENTS Core) 
include( ${CGAL_USE_FILE} )

#Inclusion of Boost
find_package(Boost 1.57.0 COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# OpenMP
#include(FindOpenMP)
#find_package(OpenMP)

# Build type
if(NOT CMAKE_BUILD_TYPE)  # Debug by default
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel"
        FORCE)
endif()

message(STATUS "Build type ${CMAKE_BUILD_TYPE}")

# Set the CFLAGS and CXXFLAGS depending on the options the user specified.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-unused-parameter")

# Debugging CFLAGS.  Turn optimizations off; turn debugging symbols on.
if(DEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O0 -Drestrict=")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++14")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++14")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -Drestrict=")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++14")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++14")
endif(DEBUG)

#if Apple OS, set the -framework Accelerate
if(APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_CXX_FLAGS} -framework Accelerate")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -framework Accelerate")
endif()

if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


#Inclusion of public headers
include_directories(include)

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

#Command file to get all the files in the src/ and include/
file(GLOB_RECURSE source_files src*/* include/*)

set(objs
 external/umat_externalM.o external/umat_externalT.o
 PROPERTIES
 EXTERNAL_OBJECT true
 GENERATED true
)

#Add the files to the lib
add_library(simcoon SHARED ${source_files})
#link against armadillo
target_link_libraries(simcoon ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} ${CGAL_LIBRARIES})

#Define lists of executables for compilation
set (All_exe_to_compile Abaqus_apply_inter solver identification L_eff Elastic_props ODF PDF)

#Compile public executable
foreach (Exe_to_compile ${All_exe_to_compile})

	#Add the solver executable
	add_executable(${Exe_to_compile} software/${Exe_to_compile}.cpp)

	#Link the solver executable with armadillo
	target_link_libraries(${Exe_to_compile} simcoon ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} ${CGAL_LIBRARIES})

endforeach (Exe_to_compile ${All_exe_to_compile})

##Testing
#Test files are in a separate source directory called test
file(GLOB_RECURSE TEST_SRCS test/*.cpp)

#Run through each source
foreach(testSrc ${TEST_SRCS})

        get_filename_component(testFolder ${testSrc} REALPATH ${CMAKE_CURRENT_SOURCE_DIR})

	file(RELATIVE_PATH testRel ${CMAKE_CURRENT_SOURCE_DIR}/test ${testFolder})
        #Extract the filename without an extension (NAME_WE)
        get_filename_component(testName ${testRel} NAME_WE)
        get_filename_component(testFolder ${testRel} DIRECTORY)

        #Add compile target
	link_directories ( ${BOOST_LIBRARY_DIRS} )
        add_executable(${testName} ${testSrc})

        #link to Boost libraries AND your targets and dependencies
        target_link_libraries(${testName} simcoon ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} ${CGAL_LIBRARIES})

        #I like to move testing binaries into a testBin directory
        set_target_properties(${testName} PROPERTIES 
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder})

        #Finally add it to test execution - 
        #Notice the WORKING_DIRECTORY and COMMAND
        add_test(NAME ${testFolder}/${testName} 
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}/${testName} )
endforeach(testSrc)


#add_library(umat OBJECT software/umat_single.cpp)
#set_property(TARGET umat PROPERTY COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -c -openm -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER")

#add_library(umatT OBJECT software/umat_singleT.cpp)
#set_property(TARGET umatT PROPERTY COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -c -openm -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER")


################################################################################
# INSTALL CONFIGURATION

if(NOT INSTALL_INCLUDE_DIR)
  set(INSTALL_INCLUDE_DIR "include")
endif()

if(NOT INSTALL_LIB_DIR)
  set(INSTALL_LIB_DIR "lib")
endif()

# executables destination
if(NOT INSTALL_BIN_DIR)
  set(INSTALL_BIN_DIR "bin")
endif()

message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
message(STATUS "INSTALL_LIB_DIR      = ${INSTALL_LIB_DIR}"     )
message(STATUS "INSTALL_INCLUDE_DIR  = ${INSTALL_INCLUDE_DIR}" )
message(STATUS "INSTALL_BIN_DIR      = ${INSTALL_BIN_DIR}"    )

install(DIRECTORY include/ DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS simcoon DESTINATION ${INSTALL_LIB_DIR})
