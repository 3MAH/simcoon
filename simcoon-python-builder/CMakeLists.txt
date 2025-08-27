# CMakeLists.txt

cmake_minimum_required(VERSION 3.15)

project(simcoon)

#set(CMAKE_VERBOSE_MAKEFILE ON)

# The version number
set(VERSION_MAJOR 1)
set(VERSION_MINOR 9)
set(VERSION_PATCH 4)

if(VERSION_PATCH MATCHES "0")
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}")
else()
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()

set(CMAKE_DISABLE_PRECOMPILE_HEADERS ON)

##### RPATH HANDLING ####
#set(CMAKE_MACOSX_RPATH 1)
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
################


# put our local cmake find scripts at the beginning of the cmake
# module search path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (MSVC) #if windows using visual c++

  # Specify the output directory for each configuration (Debug, Release, etc.)
  foreach(CONFIG_TYPE ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG_TYPE} CONFIG)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/bin")
  endforeach()
endif()

#Set executable files and library files
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

find_package(Python3 COMPONENTS Interpreter Development NumPy)

# Fetch carma if not found
include(FetchContent)
find_package(carma CONFIG QUIET)
if(NOT carma_FOUND)
  message(STATUS "carma not found, fetching from GitHub")
  FetchContent_Declare(
    carma
    GIT_REPOSITORY https://github.com/RUrlus/carma.git
    GIT_TAG        stable
  )
  FetchContent_MakeAvailable(carma)
endif()

#Inclusion of Armadillo
if (MSVC)
  find_package(pybind11 REQUIRED)
  find_package(blas REQUIRED)
  find_package(lapack REQUIRED)
  set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})  
else()
  find_package(pybind11 REQUIRED)
  find_package(Armadillo 12.6 REQUIRED)
endif()

if (NOT DEFINED USE_OPENMP)
  set(USE_OPENMP True)
endif()

if (USE_OPENMP)
  find_package(OpenMP)
endif()
find_package(simcoon REQUIRED)

#${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}

include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})
include_directories(SYSTEM ${SIMCOON_INCLUDE_DIRS})
include_directories(SYSTEM ${CARMA_INCLUDE_DIR})

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "SIMCOON_INCLUDE_DIRS = ${SIMCOON_INCLUDE_DIRS}")
message(STATUS "SIMCOON_LIBRARIES = ${SIMCOON_LIBRARIES}")
message(STATUS "ARMADILLO_INCLUDE_DIRS = ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "ARMADILLO_LIBRARIES = ${ARMADILLO_LIBRARIES}")
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
  if (USE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict= /openmp:llvm /Y-") #\Y- to disable precompile header
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict= /openmp /Y-") #\Y- to disable precompile header (don't work for an unknown reason)
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict=")
  endif()
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

#Inclusion of public headers
include_directories(include)

#Command file to get all the files in the src/ and include/
file(GLOB_RECURSE source_files_simmit src/python_wrappers/* include/simcoon/python_wrappers/*)
file(GLOB_RECURSE test_files_carma test/Carma/*.cpp test/Carma/*.hpp)

################################################################################
# SIMCOON PYBIND11 modules
#Add the files to the lib
pybind11_add_module(simmit ${source_files_simmit})

# Ensure simmit.so can always find libsimcoon.so in the same directory
if(APPLE)
  set_target_properties(simmit PROPERTIES
      INSTALL_RPATH "@loader_path"
      BUILD_WITH_INSTALL_RPATH ON
  )
elseif(UNIX)
  set_target_properties(simmit PROPERTIES
      INSTALL_RPATH "$ORIGIN"
      BUILD_WITH_INSTALL_RPATH ON
  )
endif()

if (USE_OPENMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(simmit PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

if (MSVC)
  set_target_properties(simmit PROPERTIES PREFIX "" SUFFIX ".pyd")
else()
  set_target_properties(simmit PROPERTIES PREFIX "" SUFFIX ".so")
endif()
#link against simcoon (C++ lib) and carma
target_link_libraries(simmit PRIVATE ${ARMADILLO_LIBRARIES} carma::carma Python3::Module ${SIMCOON_LIBRARIES})
set_target_properties(simmit PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

# Install the Python extension
install(TARGETS simmit
        DESTINATION ${CMAKE_SOURCE_DIR}/../python-setup/simcoon
        COMPONENT python)

# Install the core C++ library
install(TARGETS simcoon
        DESTINATION ${CMAKE_SOURCE_DIR}/../python-setup/simcoon
        COMPONENT python)


#Setup CMake to run tests
enable_testing()

if (MSVC) #if windows using visual c++

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/Carma)

  # Specify the output directory for each configuration (Debug, Release, etc.)
  foreach(CONFIG_TYPE ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG_TYPE} CONFIG)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_CURRENT_SOURCE_DIR}/test/Carma")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_CURRENT_SOURCE_DIR}/test/Carma")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_CURRENT_SOURCE_DIR}/test/Carma")
  endforeach()
endif()

set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/test)
pybind11_add_module(test_carma ${test_files_carma})

# Set output directories based on the build type
if(MSVC)
  set_target_properties(test_carma PROPERTIES
    PREFIX ""
    SUFFIX ".pyd"
  )
else()
  # For non-MSVC platforms
  set_target_properties(test_carma PROPERTIES
    PREFIX ""
    SUFFIX ".so"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test/Carma"
  )
endif()
#Link libraries
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

if(NOT INSTALL_INCLUDE_DIR)
  set(INSTALL_INCLUDE_DIR "include")
endif()

install(DIRECTORY include/ DESTINATION ${INSTALL_INCLUDE_DIR})
