# CMakeLists.txt

cmake_minimum_required(VERSION 3.15)

project(simcoon LANGUAGES C CXX)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)

# Ensure Python finds `simcoon` in the correct location
execute_process(
    COMMAND python -c "import sys; sys.path.insert(0, 'python-setup'); from simcoon.__version__ import __version__; print(__version__)"
    OUTPUT_VARIABLE PROJECT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

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
enable_testing()

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

find_package(GTest REQUIRED)

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
find_package(Boost 1.57.0 REQUIRED COMPONENTS system filesystem atomic)
message(STATUS "Boost_FOUND = ${Boost_FOUND}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
message(STATUS "Boost targets available:")
if(TARGET Boost::system)
  message(STATUS "  Boost::system")
endif()
if(TARGET Boost::filesystem) 
  message(STATUS "  Boost::filesystem")
endif()
if(TARGET Boost::atomic)
  message(STATUS "  Boost::atomic")
endif()
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

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
#link against armadillo and boost
if (MSVC AND USE_CARMA)
  target_link_libraries(simcoon ${ARMADILLO_LIBRARIES} carma::carma Boost::system Boost::filesystem Boost::atomic)
else()
  target_link_libraries(simcoon ${ARMADILLO_LIBRARIES} Boost::system Boost::filesystem Boost::atomic)
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
    add_dependencies(${testName} simcoon)
    add_custom_command(TARGET ${testName} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_BINARY_DIR}/bin
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

# Install runtime DLLs on Windows
if(WIN32)
  # Install all DLLs from build/bin to install/bin for complete installation
  install(DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
          DESTINATION ${INSTALL_BIN_DIR}
          FILES_MATCHING PATTERN "*.dll")

  # Install boost runtime dependencies by finding their actual locations
  if(Boost_FOUND)
    foreach(boost_target Boost::system Boost::filesystem Boost::atomic)
      if(TARGET ${boost_target})
        # Try to get the imported location for the current configuration
        get_target_property(boost_location ${boost_target} IMPORTED_LOCATION_RELEASE)
        if(NOT boost_location)
          get_target_property(boost_location ${boost_target} IMPORTED_LOCATION)
        endif()
        if(boost_location)
          get_filename_component(boost_dir ${boost_location} DIRECTORY)
          get_filename_component(boost_name ${boost_location} NAME_WE)
          set(boost_dll "${boost_dir}/${boost_name}.dll")
          if(EXISTS "${boost_dll}")
            install(FILES "${boost_dll}" DESTINATION ${INSTALL_BIN_DIR})
            message(STATUS "Will install boost DLL: ${boost_dll}")
          endif()
        endif()
      endif()
    endforeach()
  endif()
endif()
