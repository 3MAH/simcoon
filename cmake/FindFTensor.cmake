#.rst:
# FindFTensor
# -------------
#
# Find FTensor
#
# Find the FTensor C++ library
#
# Using FTensor:
#
# ::
#
#   find_package(FTensor REQUIRED)
#   include_directories(${FTENSOR_INCLUDE_DIRS})
#   add_executable(foo foo.cc)
#   target_link_libraries(foo ${FTENSOR_LIBRARIES})
#
# This module sets the following variables:
#
# ::
#
#   FTENSOR_FOUND - set to true if the library is found
#   FTENSOR_INCLUDE_DIRS - list of required include directories

#=============================================================================
# Copyright 2016-2017-2018 Yves Chemisky <yves.chemisky@gmail.com>
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

#Look into classical UNIX and MAC paths, what follows is for Windows users
find_path(FTENSOR_INCLUDE_DIR
  NAMES FTensor
  PATHS "$ENV{ProgramFiles}/FTensor/include"
  )
#======================

if (FTENSOR_FOUND)
  set(FTENSOR_INCLUDE_DIRS ${FTENSOR_INCLUDE_DIR})
endif ()


# Hide internal variables
mark_as_advanced(
  FTENSOR_INCLUDE_DIR
)

#======================