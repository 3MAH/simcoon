# CMakeLists.txt

cmake_minimum_required(VERSION 3.15)

project(simcoon LANGUAGES C CXX)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)

# Read project version directly from version file
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/python-setup/simcoon/__version__.py" VERSION_FILE_CONTENT)
string(REGEX MATCH "__version__ = \"([^\"]+)\"" VERSION_MATCH "${VERSION_FILE_CONTENT}")
set(PROJECT_VERSION "${CMAKE_MATCH_1}")

# Print and set the project version
message(STATUS "Project version: ${PROJECT_VERSION}")
set(PROJECT_VERSION ${PROJECT_VERSION})

# Optionally use this version in packaging/configuration
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})

message(STATUS "Configuring ${PROJECT_NAME} ${PROJECT_VERSION}")

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

# We default to debugging mode for developers.
option(DEBUG "Compile with debugging information" OFF)
option(BUILD_TESTS "Build the tests" ON)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
# Build type
if(NOT CMAKE_BUILD_TYPE)  # Debug by default
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel"
        FORCE)
endif()
if (NOT DEFINED USE_CARMA)
  set(USE_CARMA True)
endif()

message(STATUS "Build type ${CMAKE_BUILD_TYPE}")
message(STATUS "USE_CARMA          = ${USE_CARMA}"         )

# put our local cmake find scripts at the beginning of the cmake
# module search path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

#Setup CMake to run tests
if(BUILD_TESTS)
  enable_testing()
endif()

#Set executable files and library files
set(EXECUTABLE_OUTPUT_PATH bin)
set(LIBRARY_OUTPUT_PATH lib)

#Inclusion of Armadillo
if (MSVC)
  set(BLA_PKGCONFIG_BLAS True)
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Armadillo 12.6 REQUIRED)

# GTest handling - only if BUILD_TESTS is ON
if(BUILD_TESTS)
  find_package(GTest QUIET)
  
  if(NOT GTest_FOUND)
    message(STATUS "GTest not found, fetching from GitHub...")
    include(FetchContent)
    
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG        release-1.12.1
    )
    
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    
    FetchContent_MakeAvailable(googletest)
    
    # Create aliases to match find_package behavior
    add_library(GTest::gtest ALIAS gtest)
    add_library(GTest::gtest_main ALIAS gtest_main)
  endif()
endif()

if (MSVC)
  if (USE_CARMA)
    find_package(Carma REQUIRED)
  endif()
  set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})  
else()
  include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})  
endif() 

#Inclusion of Boost
find_package(Boost 1.57.0)
if (Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
endif()

# OpenMP
#include(FindOpenMP)
#find_package(OpenMP)

# Set the CFLAGS and CXXFLAGS depending on the options the user specified.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-unused-parameter")

if (MSVC) #if windows using visual c++
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict=") #\Y- to disable precompile header (don't work for an unknown reason)
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /Drestrict=")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}") #do nothing but kept anyway if required to add some options

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

  # Specify the output directory for each configuration (Debug, Release, etc.)
  foreach(CONFIG_TYPE ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG_TYPE} CONFIG)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_BINARY_DIR}/bin")
  endforeach()

elseif (UNIX AND NOT APPLE)
  #enable Fortran
  enable_language(Fortran)

  if (DEBUG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O0 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20")
  else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O3 -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20")
  endif()
else()
  #enable Fortran
  enable_language(Fortran)

  if (DEBUG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -stdlib=libc++ -O0 -fpermissive -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
  else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -stdlib=libc++ -O3 -fpermissive -Drestrict=")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++20 -stdlib=libc++")
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
if (MSVC AND USE_CARMA)
  target_link_libraries(simcoon ${ARMADILLO_LIBRARIES} carma::carma)
else()
  target_link_libraries(simcoon ${ARMADILLO_LIBRARIES})
endif()

#Define lists of executables for compilation
set (All_exe_to_compile solver identification L_eff Elastic_props ODF PDF)

#Compile public executable
foreach (Exe_to_compile ${All_exe_to_compile})

	#Add the solver executable
	add_executable(${Exe_to_compile} software/${Exe_to_compile}.cpp)

	#Link the solver executable with armadillo
  if (MSVC AND USE_CARMA)
    target_link_libraries(${Exe_to_compile} simcoon ${ARMADILLO_LIBRARIES} carma::carma)
  else()
    target_link_libraries(${Exe_to_compile} simcoon ${ARMADILLO_LIBRARIES})
  endif()
  if(NOT WIN32)
    target_link_libraries(${Exe_to_compile} -ldl)
  endif()

endforeach (Exe_to_compile ${All_exe_to_compile})

##Testing
if(BUILD_TESTS)
  #Test files are in a separate source directory called test
  file(GLOB_RECURSE TEST_SRCS test/*.cpp)
  file(GLOB_RECURSE TEST_EXTERN test_extern/*.cpp)

if(NOT MSVC)
  add_library(umat_plugin_aba SHARED testBin/Umats/UMABA/external/UMAT_ABAQUS_ELASTIC.for testBin/Umats/UMABA/external/umat_plugin_aba.cpp)
  set_target_properties(umat_plugin_aba PROPERTIES PREFIX "" SUFFIX "" LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMABA/external)
  target_link_libraries(umat_plugin_aba simcoon ${ARMADILLO_LIBRARIES})
endif()

if(MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMEXT/external)
  # Specify the output directory for each configuration (Debug, Release, etc.)
  foreach(CONFIG_TYPE ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${CONFIG_TYPE} CONFIG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} "${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMEXT/external")
  endforeach()
endif()
add_library(umat_plugin_ext SHARED testBin/Umats/UMEXT/external/umat_plugin_ext.cpp)
if(MSVC)
  set_target_properties(umat_plugin_ext PROPERTIES PREFIX "" SUFFIX ".dll")
  if (USE_CARMA)
    target_link_libraries(umat_plugin_ext simcoon ${ARMADILLO_LIBRARIES} carma::carma)
  else()
    target_link_libraries(umat_plugin_ext simcoon ${ARMADILLO_LIBRARIES})
  endif()
else()
  set_target_properties(umat_plugin_ext PROPERTIES PREFIX "" SUFFIX "" LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/Umats/UMEXT/external)
  target_link_libraries(umat_plugin_ext simcoon ${ARMADILLO_LIBRARIES})
endif()

#Run through each source
foreach(testSrc ${TEST_SRCS})

  get_filename_component(testFolder ${testSrc} REALPATH ${CMAKE_CURRENT_SOURCE_DIR})
  file(RELATIVE_PATH testRel ${CMAKE_CURRENT_SOURCE_DIR}/test ${testFolder})

  # Extract the filename without an extension (NAME_WE)
  get_filename_component(testName ${testRel} NAME_WE)
  get_filename_component(testFolder ${testRel} DIRECTORY)

  #Add compile target
  add_executable(${testName} ${testSrc})

  if (MSVC AND USE_CARMA)
  	target_link_libraries(${testName} simcoon ${ARMADILLO_LIBRARIES} GTest::gtest GTest::gtest_main carma::carma)
  else()
  	target_link_libraries(${testName} simcoon ${ARMADILLO_LIBRARIES} GTest::gtest GTest::gtest_main)
  endif()
  target_include_directories(${testName} PRIVATE ${ARMADILLO_INCLUDE_DIRS} ${GTEST_INCLUDE_DIRS})    

  # Custom runtime output directory for test binaries
  set_target_properties(${testName} PROPERTIES 
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
      RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
      RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder})

  if(MSVC)
    add_custom_command(TARGET ${testName} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_BINARY_DIR}/bin/simcoon.dll  
        ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder})
  endif()

  # Add test with proper working directory and command
  add_test(NAME ${testFolder}/${testName} 
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}/${testName})        
endforeach(testSrc)

if(NOT MSVC)
  #Run through each source
  foreach(testSrc ${TEST_EXTERN})

    get_filename_component(testFolder ${testSrc} REALPATH ${CMAKE_CURRENT_SOURCE_DIR})
    file(RELATIVE_PATH testRel ${CMAKE_CURRENT_SOURCE_DIR}/test_extern ${testFolder})

    # Extract the filename without an extension (NAME_WE)
    get_filename_component(testName ${testRel} NAME_WE)
    get_filename_component(testFolder ${testRel} DIRECTORY)

    #Add compile target
    add_executable(${testName} ${testSrc})

    target_link_libraries(${testName} simcoon GTest::gtest GTest::gtest_main)
    target_include_directories(${testName} PRIVATE ${Boost_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ${GTEST_INCLUDE_DIRS})        

    # Custom runtime output directory for test binaries
    set_target_properties(${testName} PROPERTIES 
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder})

    # Add test with proper working directory and command
    add_test(NAME ${testFolder}/${testName} 
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testFolder}/${testName})        
  endforeach(testSrc)
endif()

endif() # BUILD_TESTS

#set_property(TARGET FOR PROPERTY)
#add_library(ABA OBJECT testBin/Umats/UMABA/external/umat_plugin_aba.cpp)
#set_property(TARGET ABA PROPERTY COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fPIC o umat_plugin_aba -lsimcoon -larmadillo -lgfortran")

#add_library(umat OBJECT software/umat_single.cpp)
#set_property(TARGET umat PROPERTY COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -c -openm -larmadillo -std=c++20 -DARMA_DONT_USE_WRAPPER")

#add_library(umatT OBJECT software/umat_singleT.cpp)
#set_property(TARGET umatT PROPERTY COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -c -openm -larmadillo -std=c++20 -DARMA_DONT_USE_WRAPPER")


################################################################################
# INSTALL CONFIGURATION

set(INSTALL_INCLUDE_DIR include)
set(INSTALL_LIB_DIR lib)

# executables destination
if(NOT INSTALL_BIN_DIR)
  set(INSTALL_BIN_DIR bin)
endif()

message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
message(STATUS "INSTALL_LIB_DIR      = ${INSTALL_LIB_DIR}"     )
message(STATUS "INSTALL_INCLUDE_DIR  = ${INSTALL_INCLUDE_DIR}" )
message(STATUS "INSTALL_BIN_DIR      = ${INSTALL_BIN_DIR}"    )

install(DIRECTORY include/ DESTINATION ${INSTALL_INCLUDE_DIR})
install(TARGETS simcoon)

################################################################################
# PYTHON BINDINGS CONFIGURATION

if(BUILD_PYTHON_BINDINGS)
  message(STATUS "Building Python bindings...")
  
  # Find Python dependencies for bindings
  find_package(Python3 COMPONENTS Interpreter Development QUIET)
  if(NOT Python3_FOUND)
    message(WARNING "Python3 not found. Python bindings will be disabled.")
    set(BUILD_PYTHON_BINDINGS OFF)
  else()
    message(STATUS "Python3 found: ${Python3_EXECUTABLE}")
    
    # Try to find pybind11, first with standard search
    find_package(pybind11 QUIET)
    if(NOT pybind11_FOUND)
      message(STATUS "pybind11 not found in standard locations, trying to locate via Python...")
      
      # Try to get pybind11 cmake directory from Python
      execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
        OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE PYBIND11_PYTHON_RESULT
      )
      
      if(PYBIND11_PYTHON_RESULT EQUAL 0 AND PYBIND11_CMAKE_DIR)
        message(STATUS "Found pybind11 cmake directory via Python: ${PYBIND11_CMAKE_DIR}")
        
        # Add to CMAKE_PREFIX_PATH and try again
        list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_CMAKE_DIR}")
        find_package(pybind11 QUIET)
        
        if(pybind11_FOUND)
          message(STATUS "pybind11 found via Python detection")
        else()
          message(WARNING "pybind11 cmake directory found but package not loaded correctly")
          set(BUILD_PYTHON_BINDINGS OFF)
        endif()
      else()
        message(WARNING "pybind11 not found. Python bindings require pybind11.")
        message(STATUS "Please install: pip install pybind11")
        set(BUILD_PYTHON_BINDINGS OFF)
      endif()
    else()
      message(STATUS "pybind11 found")
    endif()
    
    # If pybind11 was found (either way), set up the rest
    if(pybind11_FOUND)
      # Determine pip command based on Python environment
      string(FIND "${Python3_EXECUTABLE}" ".venv" UV_DETECTED)
      if(UV_DETECTED GREATER -1)
        set(PIP_COMMAND "uv" "pip")
        message(STATUS "Detected uv environment, using 'uv pip'")
      else()
        set(PIP_COMMAND "${Python3_EXECUTABLE}" "-m" "pip")
        message(STATUS "Using standard pip")
      endif()
      
      # Try to find carma, if not found, clone it
      find_package(carma CONFIG QUIET)
      if(NOT carma_FOUND)
        message(STATUS "carma not found, fetching from GitHub...")
        include(FetchContent)
        
        FetchContent_Declare(
          carma
          GIT_REPOSITORY https://github.com/RUrlus/carma.git
          GIT_TAG        v0.8.0
        )
        
        FetchContent_MakeAvailable(carma)
        message(STATUS "carma fetched and configured")
      else()
        message(STATUS "carma found")
      endif()
      
      # Set variables for simcoon-python-builder
      set(SIMCOON_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
      set(SIMCOON_LIBRARIES simcoon)
      
      # Add simcoon-python-builder as subdirectory
      add_subdirectory(simcoon-python-builder)
      
      # Set RPATH for simmit to find libsimcoon.dylib in the same directory
      if(APPLE)
        set_target_properties(simmit PROPERTIES
          INSTALL_RPATH "@loader_path"
          BUILD_WITH_INSTALL_RPATH TRUE)
      else()
        set_target_properties(simmit PROPERTIES
          INSTALL_RPATH "\$ORIGIN"
          BUILD_WITH_INSTALL_RPATH TRUE)
      endif()
    endif()
  endif()
endif()

if(BUILD_PYTHON_BINDINGS)
  # Install Python module and required libraries
  install(TARGETS simmit simcoon DESTINATION simcoon COMPONENT python)
  
  # Install Python package files
  install(DIRECTORY python-setup/simcoon/ 
          DESTINATION simcoon 
          COMPONENT python
          FILES_MATCHING PATTERN "*.py")
endif()
