# CMakeLists.txt

cmake_minimum_required(VERSION 3.15)
project(simcoon CXX C)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)

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

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

# We default to debugging mode for developers.
option(DEBUG "Compile with debugging information" OFF)
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

#Set executable files and library files
if (MSVC)
  set(EXECUTABLE_OUTPUT_PATH bin)
  set(LIBRARY_OUTPUT_PATH lib)
else()
  set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
  set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})
endif()

#Inclusion of Armadillo
if (MSVC)
  find_package(pybind11 REQUIRED)
  find_package(carma CONFIG REQUIRED)
  find_package(blas REQUIRED)
  find_package(lapack REQUIRED)
  set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})    
  include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
else()
  include(FindArmadillo)
  find_package(Armadillo 5.2 REQUIRED)
  include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
endif()

#Find CGAL
find_package(CGAL COMPONENTS Core REQUIRED)
#include_directories(SYSTEM ${CGAL_INCLUDE_DIRS})

#Inclusion of Boost
find_package(Boost 1.57.0 COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# OpenMP
#include(FindOpenMP)
#find_package(OpenMP)

# Set the CFLAGS and CXXFLAGS depending on the options the user specified.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-unused-parameter")

if (MSVC) #if windows using visual c++
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict= /Y-") #\Y- to disable precompile header (don't work for an unknown reason)
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict=")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
elseif (UNIX AND NOT APPLE)
  #enable Fortran
  enable_language(Fortran)

  if (DEBUG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O0 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++14")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++14")
  else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++14")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++14")
  endif()
else()
  #enable Fortran
  enable_language(Fortran)

  if (DEBUG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -stdlib=libc++ -O0 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++14 -stdlib=libc++")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++14 -stdlib=libc++")
  else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -stdlib=libc++ -O3 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++14 -stdlib=libc++")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++14 -stdlib=libc++")
  endif()
endif()

if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

message(STATUS "CMAKE_SYSTEM_NAME          = ${CMAKE_SYSTEM_NAME}"         )
message(STATUS "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}"     )
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_CXX_FLAGS            = ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_COMPILER_IS_GNUCXX   = ${CMAKE_COMPILER_IS_GNUCXX}"  )
message(STATUS "BUILD_SHARED_LIBS          = ${BUILD_SHARED_LIBS}"         )
message(STATUS "DETECT_HDF5                = ${DETECT_HDF5}"               )
message(STATUS "ARMADILLO_INCLUDE_DIRS = ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "ARMADILLO_LIBRARIES = ${ARMADILLO_LIBRARIES}")


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
if (MSVC)
  #for windows, need to include carma to avoid bug in the pybhon-builder (https://github.com/RUrlus/carma/issues/91)
  target_link_libraries(simcoon ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} carma::carma CGAL::CGAL CGAL::CGAL_Core)
else()
  target_link_libraries(simcoon ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core)
endif()

#Define lists of executables for compilation
set (All_exe_to_compile Abaqus_apply_inter Salome_apply_inter Generic_apply_inter  extract_Abaqus solver identification L_eff Elastic_props ODF PDF)

#Compile public executable
foreach (Exe_to_compile ${All_exe_to_compile})

	#Add the solver executable
	add_executable(${Exe_to_compile} software/${Exe_to_compile}.cpp)

	#Link the solver executable with armadillo
	target_link_libraries(${Exe_to_compile} simcoon ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core -ldl)

endforeach (Exe_to_compile ${All_exe_to_compile})

##Testing
#Test files are in a separate source directory called test
file(GLOB_RECURSE TEST_SRCS test/*.cpp)

add_library(umat_plugin_aba SHARED testBin/Umats/UMABA/external/UMAT_ABAQUS_ELASTIC.for testBin/Umats/UMABA/external/umat_plugin_aba.cpp)
set_target_properties(umat_plugin_aba PROPERTIES PREFIX "" SUFFIX "" LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMABA/external)
target_link_libraries(umat_plugin_aba simcoon  ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core)

add_library(umat_plugin_ext SHARED testBin/Umats/UMEXT/external/umat_plugin_ext.cpp)
set_target_properties(umat_plugin_ext PROPERTIES PREFIX "" SUFFIX "" LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMEXT/external)
target_link_libraries(umat_plugin_ext simcoon  ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core)

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
        target_link_libraries(${testName} simcoon ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core -ldl)

        #I like to move testing binaries into a testBin directory
        set_target_properties(${testName} PROPERTIES 
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder})

        #Finally add it to test execution - 
        #Notice the WORKING_DIRECTORY and COMMAND
        add_test(NAME ${testFolder}/${testName} 
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}/${testName} )
endforeach(testSrc)

#set_property(TARGET FOR PROPERTY)
#add_library(ABA OBJECT testBin/Umats/UMABA/external/umat_plugin_aba.cpp)
#set_property(TARGET ABA PROPERTY COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fPIC o umat_plugin_aba -lsimcoon -larmadillo -lgfortran")

#add_library(umat OBJECT software/umat_single.cpp)
#set_property(TARGET umat PROPERTY COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -c -openm -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER")

#add_library(umatT OBJECT software/umat_singleT.cpp)
#set_property(TARGET umatT PROPERTY COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -c -openm -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER")


################################################################################
# INSTALL CONFIGURATION

set(INSTALL_INCLUDE_DIR include)
set(INSTALL_LIB_DIR lib)

# executables destination
#if(NOT INSTALL_BIN_DIR)
#  set(INSTALL_BIN_DIR "bin")
#endif()

message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
message(STATUS "INSTALL_LIB_DIR      = ${INSTALL_LIB_DIR}"     )
message(STATUS "INSTALL_INCLUDE_DIR  = ${INSTALL_INCLUDE_DIR}" )
message(STATUS "INSTALL_BIN_DIR      = ${INSTALL_BIN_DIR}"    )

install(DIRECTORY include/ DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS simcoon DESTINATION ${INSTALL_LIB_DIR})
