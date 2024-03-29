# CMakeLists.txt

cmake_minimum_required(VERSION 3.15)

project(simcoon)

#set(CMAKE_VERBOSE_MAKEFILE ON)

# The version number
set(VERSION_MAJOR 0)
set(VERSION_MINOR 9)
set(VERSION_PATCH 0)

if(VERSION_PATCH MATCHES "0")
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}")
else()
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()

##### RPATH HANDLING ####
#set(CMAKE_MACOSX_RPATH 1)
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
################


# put our local cmake find scripts at the beginning of the cmake
# module search path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib 
    CACHE PATH "Output directory for static libraries.")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    CACHE PATH "Output directory for shared libraries.")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    CACHE PATH "Output directory for executables and DLL's.")

#Set executable files and library files
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

find_package(Python3 COMPONENTS Interpreter Development NumPy)

if (Python3_FOUND)
  if(MSVC)
    set(CMAKE_FIND_FRAMEWORK NEVER)
      FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX} numpy${PYTHON_VERSION_SUFFIX} REQUIRED)
      #FIND_PACKAGE(PythonInterp 3)
      #FIND_PACKAGE(PythonLibs 3 REQUIRED)
      #FIND_PACKAGE(NumPy REQUIRED)
  elseif (UNIX AND NOT APPLE)
    FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX} numpy${PYTHON_VERSION_SUFFIX} REQUIRED)
    FIND_PACKAGE(PythonInterp 3)
    FIND_PACKAGE(PythonLibs 3 REQUIRED)
    FIND_PACKAGE(NumPy REQUIRED)
  else()
    set(CMAKE_FIND_FRAMEWORK NEVER)
      FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} numpy${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED)
      FIND_PACKAGE(PythonInterp 3)
      FIND_PACKAGE(PythonLibs 3 REQUIRED)
      FIND_PACKAGE(NumPy REQUIRED)
  endif()
else()
    message("Python not found")
endif()

#Inclusion of Armadillo
if (MSVC)
  find_package(pybind11 REQUIRED)
  find_package(carma CONFIG REQUIRED)
  find_package(blas REQUIRED)
  find_package(lapack REQUIRED)
  set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})  
else()
  find_package(pybind11 REQUIRED)
  find_package(Armadillo 12.6 REQUIRED)
  find_package(carma CONFIG REQUIRED)
endif()

find_package(simcoon REQUIRED)
find_package(arma2numpy REQUIRED)

#${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}

include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
include_directories(SYSTEM ${Boost_INCLUDE_DIRS} )
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})
include_directories(SYSTEM ${SIMCOON_INCLUDE_DIRS})
include_directories(SYSTEM ${CARMA_INCLUDE_DIR})

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
message(STATUS "SIMCOON_INCLUDE_DIRS = ${SIMCOON_INCLUDE_DIRS}")
message(STATUS "SIMCOON_LIBRARIES = ${SIMCOON_LIBRARIES}")
message(STATUS "ARMADILLO_INCLUDE_DIRS = ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "ARMADILLO_LIBRARIES = ${ARMADILLO_LIBRARIES}")
message(STATUS "arma2numpy_LIBRARIES = ${arma2numpy_LIBRARIES}")
message(STATUS "carma_INCLUDE_DIR = ${carma_INCLUDE_DIR}")

# Build type
if(NOT CMAKE_BUILD_TYPE)  # Debug by default
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel"
        FORCE)
endif()

# Set the CFLAGS and CXXFLAGS depending on the options the user specified.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-unused-parameter")


if (MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict= /Y-") #\Y- to disable precompile header (don't work for an unknown reason)
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict=")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
elseif (UNIX AND NOT APPLE)
  if (DEBUG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O0 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++20")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20")
  else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O3 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++20")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20")
  endif()
else()
  if (DEBUG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -stdlib=libc++ -O0 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
  else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -stdlib=libc++ -O3 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
  endif()
endif()

if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#Inclusion of public headers
include_directories(include)

#Command file to get all the files in the src/ and include/
file(GLOB_RECURSE source_files_simmit src/python_wrappers/* include/simcoon/python_wrappers/*)
file(GLOB_RECURSE test_files_carma test/Carma/*.cpp test/Carma/*.hpp)

#Setup CMake to run tests
enable_testing()

#Test Carma, link against carma, 
Python3_add_library(test_carma MODULE ${test_files_carma})
set(LIBRARY_OUTPUT_PATH {CMAKE_CURRENT_SOURCE_DIR}/test)
set_target_properties(test_carma PROPERTIES PREFIX "" SUFFIX ".so" LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma)
target_link_libraries(test_carma PRIVATE carma::carma ${ARMADILLO_LIBRARIES} Python3::Module)

#Finally add it to test execution -
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME test_carma_roundtrip WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma
		COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma/run_test_roundtrip.py)
add_test(NAME test_carma_arr_to_mat WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma
		COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma/run_test_arr_to_mat.py)    

#Finally add it to test execution -
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME simmit_test WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/simmit_test
		COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/simmit_test/run_test.py)

################################################################################
# SIMCOON modules
#Add the files to the lib
#add_library(simmitboost MODULE ${source_files_simmit})
Python3_add_library(simmitboost MODULE ${source_files_simmit})
#Wrapper library set_target properties
if (MSVC)
  set_target_properties(simmitboost PROPERTIES PREFIX "" SUFFIX ".pyd")
else()
  set_target_properties(simmitboost PROPERTIES PREFIX "" SUFFIX ".so")
endif()
#link against simcoon (C++ lib) and arma2numpy, and Boost for the wrappers
target_link_libraries(simmitboost PRIVATE ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} carma::carma Python3::Module ${SIMCOON_LIBRARIES} ${arma2numpy_LIBRARIES} )
#if (UNIX AND APPLE)
#	set_target_properties(simmitboost PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
#endif()

################################################################################

################################################################################
# SIMCOON PYBIND11 modules
#Add the files to the lib
pybind11_add_module(simmit ${source_files_simmit})
#Wrapper library set_target properties
if (MSVC)
  set_target_properties(simmit PROPERTIES PREFIX "" SUFFIX ".pyd")
else()
  set_target_properties(simmit PROPERTIES PREFIX "" SUFFIX ".so")
endif()
#link against simcoon (C++ lib) and arma2numpy, and Boost for the wrappers
target_link_libraries(simmit PRIVATE ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} carma::carma Python3::Module ${SIMCOON_LIBRARIES} ${arma2numpy_LIBRARIES})
#if (UNIX AND APPLE)
#	set_target_properties(simmit PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
#endif()

################################################################################

if(NOT INSTALL_INCLUDE_DIR)
  set(INSTALL_INCLUDE_DIR "include")
endif()

install(DIRECTORY include/ DESTINATION ${INSTALL_INCLUDE_DIR})
