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
#   include_directories(${ARMA2NUMPY_INCLUDE_DIRS})
#   add_library(foo MODULE ${source_files})
#   target_link_libraries(foo ${ARMA2NUMPY_LIBRARIES})
#
# This module sets the following variables:
#
# ::
#
#   ARMA2NUMPY_FOUND - set to true if the library is found
#   ARMA2NUMPY_INCLUDE_DIRS - list of required include directories
#   ARMA2NUMPY_LIBRARIES - list of libraries to be linked

#=============================================================================
# Copyright 2016 Yves Chemisky <yves.chemisky@gmail.com>
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

#Look into classical UNIX and MAC paths, what follows is for Windows users
find_library(ARMA2NUMPY_LIBRARY
  NAMES arma2numpy
  PATHS "$ENV{ProgramFiles}/arma2numpy/lib"  "$ENV{ProgramFiles}/arma2numpy/lib64" "$ENV{ProgramFiles}/arma2numpy"
  )
find_path(ARMA2NUMPY_INCLUDE_DIR
  NAMES arma2numpy
  PATHS "$ENV{ProgramFiles}/arma2numpy/include"
  )
#======================

if (ARMA2NUMPY_FOUND)
  set(ARMA2NUMPY_INCLUDE_DIRS ${ARMA2NUMPY_INCLUDE_DIR})
  set(ARMA2NUMPY_LIBRARIES ${ARMA2NUMPY_LIBRARY})
endif ()

# Hide internal variables
mark_as_advanced(
  SMARTPLUS_INCLUDE_DIR
  SMARTPLUS_LIBRARY)

#======================