#.rst:
# Findsimcoon
# -------------
#
# Find simcoon
#
# Find the simcoon C++ library
#
# Using simcoon:
#
# ::
#
#   find_package(simcoon REQUIRED)
#   include_directories(${SIMCOON_INCLUDE_DIRS})
#   add_executable(foo foo.cc)
#   target_link_libraries(foo ${SIMCOON_LIBRARIES})
#
# This module sets the following variables:
#
# ::
#
#   SIMCOON_FOUND - set to true if the library is found
#   SIMCOON_INCLUDE_DIRS - list of required include directories
#   SIMCOON_LIBRARIES - list of libraries to be linked

#=============================================================================
# Copyright 2016-2017-2018 Yves Chemisky <yves.chemisky@gmail.com>
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

# UNIX paths are standard, no need to write.
find_library(SIMCOON_LIBRARY
  NAMES simcoon
  PATHS "$ENV{ProgramFiles}/simcoon/lib"  "$ENV{ProgramFiles}/simcoon/lib64" "$ENV{ProgramFiles}/simcoon"
  )
find_path(SIMCOON_INCLUDE_DIR
  NAMES simcoon
  PATHS "$ENV{ProgramFiles}/simcoon/include"
  )

include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
find_package_handle_standard_args(simcoon
  REQUIRED_VARS SIMCOON_LIBRARY SIMCOON_INCLUDE_DIR)
#  VERSION_VAR SIMCOON_VERSION_STRING)
# version_var fails with cmake < 2.8.4.

if (SIMCOON_FOUND)
  set(SIMCOON_INCLUDE_DIRS ${SIMCOON_INCLUDE_DIR})
  set(SIMCOON_LIBRARIES ${SIMCOON_LIBRARY})
endif ()

# Hide internal variables
mark_as_advanced(
  SIMCOON_INCLUDE_DIR
  SIMCOON_LIBRARY)
