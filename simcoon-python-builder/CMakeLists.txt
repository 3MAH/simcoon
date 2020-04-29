# CMakeLists.txt

cmake_minimum_required(VERSION 2.6)

project(simcoon)

# The version number
set(VERSION_MAJOR 0)
set(VERSION_MINOR 9)
set(VERSION_PATCH 0)

if(VERSION_PATCH MATCHES "0")
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}")
else()
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()

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

FIND_PACKAGE(PythonInterp)

if (PYTHONINTERP_FOUND)
  if (UNIX AND NOT APPLE)
    if (PYTHON_VERSION_MAJOR EQUAL 3)
        FIND_PACKAGE(CGAL REQUIRED COMPONENTS Core)
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX} numpy${PYTHON_VERSION_SUFFIX} REQUIRED)
        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3 REQUIRED)
        FIND_PACKAGE(NumPy REQUIRED)
    else()
        FIND_PACKAGE(CGAL REQUIRED COMPONENTS Core)
        FIND_PACKAGE(Boost COMPONENTS python numpy REQUIRED)
        FIND_PACKAGE(PythonInterp)
        FIND_PACKAGE(PythonLibs REQUIRED)
        FIND_PACKAGE(NumPy REQUIRED)
    endif()
  else()
    set(CMAKE_FIND_FRAMEWORK NEVER)
    if (PYTHON_VERSION_MAJOR EQUAL 3)
        FIND_PACKAGE(CGAL REQUIRED COMPONENTS Core)
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} numpy${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED)
        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3 REQUIRED)
        FIND_PACKAGE(NumPy REQUIRED)
    else()
        FIND_PACKAGE(CGAL REQUIRED COMPONENTS Core)
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} numpy${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED)
        FIND_PACKAGE(PythonInterp)
        FIND_PACKAGE(PythonLibs REQUIRED)
        FIND_PACKAGE(NumPy REQUIRED)
    endif()
  endif()
else()
    message("Python not found")
endif()

#${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}

include_directories(SYSTEM ${Boost_INCLUDE_DIRS} )
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

message(STATUS "BOOST_INCLUDE_DIRS = ${BOOST_INCLUDE_DIRS}")
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")

find_package(simcoon) 
if(SIMCOON_FOUND)
	include_directories(SYSTEM ${SIMCOON_INCLUDE_DIRS})
else()	
	FIND_LIBRARY(SIMCOON_LIBRARY NAMES  simcoon PATHS ${CMAKE_CURRENT_SOURCE_DIR}/lib/)
endif(SIMCOON_FOUND)

# Build type
if(NOT CMAKE_BUILD_TYPE)  # Debug by default
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel"
        FORCE)
endif()

# Set the CFLAGS and CXXFLAGS depending on the options the user specified.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-unused-parameter")


if (UNIX AND NOT APPLE)
  if (DEBUG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O0 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++14")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++14")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++14")
  else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++14")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++14")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++14")
  endif()
else()
  if (DEBUG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -stdlib=libc++ -O0 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++14 -stdlib=libc++")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++14 -stdlib=libc++")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++14 -stdlib=libc++")
  else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -stdlib=libc++ -O3 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++14 -stdlib=libc++")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -std=c++14 -stdlib=libc++")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++14 -stdlib=libc++")
  endif()
endif()

if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#Inclusion of public headers
include_directories(include)

#Command file to get all the files in the src/ and include/
file(GLOB_RECURSE source_files_smart src/python_wrappers/* include/simcoon/python_wrappers/*)

#Command file to get all the files in the src/ and include/
file(GLOB_RECURSE source_files_arma2numpy src/arma2numpy/* include/arma2numpy/*)

#Add the files to the lib
add_library(arma2numpy SHARED ${source_files_arma2numpy})

#link against simcoon and Armadillo, and Boost for the numpy
target_link_libraries(arma2numpy ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core ${PYTHON_LIBRARIES} ${SIMCOON_LIBRARY})

#Setup CMake to run tests
enable_testing()

##Testing
#Test files are in a separate source directory called test
file(GLOB_RECURSE test_files test/*.cpp)

#Add the files to the lib
#add_library(Tarma2numpy MODULE ${test_files})
PYTHON_ADD_MODULE(Tarma2numpy ${test_files})
#Move testing binaries into a testBin directory
#Set executable files and library files
set(LIBRARY_OUTPUT_PATH {CMAKE_CURRENT_SOURCE_DIR}/test)

set_target_properties(Tarma2numpy PROPERTIES PREFIX "" SUFFIX ".so" LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/test/arma2numpy)

#link against simcoon and arma2numpy, and Boost for the wrappers
target_link_libraries(Tarma2numpy ${Boost_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core ${PYTHON_LIBRARIES} ${ARMADILLO_LIBRARIES} arma2numpy)

#Finally add it to test execution -
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME Tarma2numpy WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/arma2numpy
		COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/arma2numpy/run_test.py)

################################################################################
# SIMCOON modules
#Add the files to the lib
add_library(simmit MODULE ${source_files_smart})
#Wrapper library set_target properties
set_target_properties(simmit PROPERTIES PREFIX "" SUFFIX ".so")
#link against simcoon (C++ lib) and arma2numpy, and Boost for the wrappers
target_link_libraries(simmit ${Boost_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core ${PYTHON_LIBRARIES} ${SIMCOON_LIBRARY} arma2numpy)


################################################################################
# INSTALL CONFIGURATION OF ARMA2NUMPY

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
install(TARGETS arma2numpy DESTINATION ${INSTALL_LIB_DIR})

