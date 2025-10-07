#.rst:
# Findarma2numpy
# -------------
#
# Find arma2numpy
#
# Find the arma2numpy C++ library
#
# Using arma2numpy:
#
# ::
#
#   find_package(arma2numpy REQUIRED)
#   include_directories(${arma2numpy_INCLUDE_DIRS})
#   add_executable(foo foo.cc)
#   target_link_libraries(foo ${arma2numpy_LIBRARIES})
#
# This module sets the following variables:
#
# ::
#
#   arma2numpy_FOUND - set to true if the library is found
#   arma2numpy_INCLUDE_DIRS - list of required include directories
#   arma2numpy_LIBRARIES - list of libraries to be linked

#=============================================================================
# Copyright 2016-2017-2018 Yves Chemisky <yves.chemisky@gmail.com>
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

# UNIX paths are standard, no need to write.
find_library(arma2numpy_LIBRARY
  NAMES arma2numpy
#  PATHS "$ENV{ProgramFiles}/arma2numpy/lib"  "$ENV{ProgramFiles}/arma2numpy/lib64" "$ENV{ProgramFiles}/arma2numpy"
  )
find_path(arma2numpy_INCLUDE_DIR
  NAMES simcoon
#  PATHS "$ENV{ProgramFiles}/simcoon/include"
  )

include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
find_package_handle_standard_args(arma2numpy
  REQUIRED_VARS arma2numpy_LIBRARY arma2numpy_INCLUDE_DIR)
#  VERSION_VAR arma2numpy_VERSION_STRING)
# version_var fails with cmake < 2.8.4.

if (arma2numpy_FOUND)
  set(arma2numpy_INCLUDE_DIRS ${arma2numpy_INCLUDE_DIR})
  set(arma2numpy_LIBRARIES ${arma2numpy_LIBRARY})
endif ()

# Hide internal variables
mark_as_advanced(
  arma2numpy_INCLUDE_DIR
  arma2numpy_LIBRARY)
