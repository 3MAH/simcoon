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
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#set(CMAKE_SKIP_BUILD_RPATH FALSE)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
################

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

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
    FIND_PACKAGE(CGAL REQUIRED COMPONENTS Core)
    FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX} numpy${PYTHON_VERSION_SUFFIX} REQUIRED)
    #FIND_PACKAGE(PythonInterp 3)
    #FIND_PACKAGE(PythonLibs 3 REQUIRED)
    #FIND_PACKAGE(NumPy REQUIRED)
  elseif (UNIX AND NOT APPLE)
    FIND_PACKAGE(CGAL REQUIRED COMPONENTS Core)
    FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX} numpy${PYTHON_VERSION_SUFFIX} REQUIRED)
    FIND_PACKAGE(PythonInterp 3)
    FIND_PACKAGE(PythonLibs 3 REQUIRED)
    FIND_PACKAGE(NumPy REQUIRED)
  else()
    set(CMAKE_FIND_FRAMEWORK NEVER)
      FIND_PACKAGE(CGAL REQUIRED COMPONENTS Core)
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
  find_package(carma CONFIG REQUIRED)
  find_package(blas REQUIRED)
  find_package(lapack REQUIRED)
  set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}) 
  include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
else()
  find_package(Armadillo 12.6 REQUIRED)
  include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
endif()

find_package(simcoon REQUIRED)

#${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}

include_directories(SYSTEM ${Boost_INCLUDE_DIRS} )
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})
include_directories(SYSTEM ${SIMCOON_INCLUDE_DIRS})
include_directories(SYSTEM ${CGAL_INCLUDE_DIRS})

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
message(STATUS "SIMCOON_INCLUDE_DIRS = ${SIMCOON_INCLUDE_DIRS}")
message(STATUS "SIMCOON_LIBRARIES = ${SIMCOON_LIBRARIES}")

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
  
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict=")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
elseif (UNIX AND NOT APPLE)
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
file(GLOB_RECURSE source_files_arma2numpy src/arma2numpy/* include/arma2numpy/*)

#Add the files to the lib
add_library(arma2numpy SHARED ${source_files_arma2numpy})

#link against simcoon and Armadillo, and Boost for the numpy
target_link_libraries(arma2numpy PRIVATE ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core ${SIMCOON_LIBRARIES} Python3::Module)
#if (UNIX AND APPLE)
#    set_target_properties(arma2numpy PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
#endif()

#Setup CMake to run tests
enable_testing()

##Testing
#Test files are in a separate source directory called test
file(GLOB_RECURSE test_files_arma2numpy test/arma2numpy/*.cpp test/arma2numpy/*.hpp)
file(GLOB_RECURSE test_files_CM_func test/CM_func/*.cpp test/CM_func/*.hpp)

#Add the files to the lib
#add_library(Tarma2numpy MODULE ${test_files_arma2numpy})
#add_library(TCM_func MODULE ${test_files_CM_func})
Python3_add_library(Tarma2numpy MODULE ${test_files_arma2numpy})
Python3_add_library(TCM_func MODULE ${test_files_CM_func})
#Move testing binaries into a testBin directory
#Set executable files and library files
set(LIBRARY_OUTPUT_PATH {CMAKE_CURRENT_SOURCE_DIR}/test)

set_target_properties(Tarma2numpy PROPERTIES PREFIX "" SUFFIX ".so" LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/test/arma2numpy)
set_target_properties(TCM_func PROPERTIES PREFIX "" SUFFIX ".so" LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/test/CM_func)

#link against simcoon and arma2numpy, and Boost for the wrappers
target_link_libraries(Tarma2numpy PRIVATE ${Boost_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core ${ARMADILLO_LIBRARIES} ${SIMCOON_LIBRARIES} Python3::Module arma2numpy)

#Finally add it to test execution -
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME Tarma2numpy WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/arma2numpy
		COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/arma2numpy/run_test.py)

#link against simcoon and arma2numpy, and Boost for the wrappers
target_link_libraries(TCM_func PRIVATE ${Boost_LIBRARIES} CGAL::CGAL CGAL::CGAL_Core ${ARMADILLO_LIBRARIES} ${SIMCOON_LIBRARIES} Python3::Module arma2numpy)

#Finally add it to test execution -
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME TCM_func WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/CM_func
		COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/CM_func/run_test.py)

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
message(STATUS "CMAKE_INSTALL_RPATH      = ${CMAKE_INSTALL_RPATH}"    )


install(DIRECTORY include/ DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS arma2numpy DESTINATION ${INSTALL_LIB_DIR})

