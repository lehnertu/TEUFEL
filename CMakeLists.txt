cmake_minimum_required (VERSION 2.6)
project (TEUFEL)

set(CMAKE_VERBOSE_MAKEFILE ON)

# Set default build type.
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "CMAKE_BUILD_TYPE not given; setting to 'Release'.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# compiler flags
add_definitions(-std=c++11 -Wall)

find_package(HDF5 COMPONENTS CXX REQUIRED)
if (HDF5_FOUND)
    message(STATUS "HDF5 includes ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5 libraries ${HDF5_LIBRARIES}")
endif (HDF5_FOUND)

find_package(ZLIB REQUIRED)

# our local SDDS library
set(SDDS_DIR ${PROJECT_SOURCE_DIR}/lib/SDDSToolKit-devel-3.4)
find_library(SDDS_LIB NAMES SDDS1 HINTS ${SDDS_DIR})
find_library(SDDS_MDBLIB NAMES mdblib HINTS ${SDDS_DIR})
find_library(SDDS_MDBCOMMON NAMES mdbcommon HINTS ${SDDS_DIR})
find_library(SDDS_LZMA NAMES lzma HINTS ${SDDS_DIR})

# where to find header files
include_directories( ${PROJECT_BINARY_DIR}
                     ${PROJECT_SOURCE_DIR}/src 
                     ${HDF5_INCLUDE_DIRS}
                     ${SDDS_DIR})

# we put all objects into a library that can be linked to the various executables
aux_source_directory(${PROJECT_SOURCE_DIR}/src lib_sources )
# the main program is not part of the library
list(REMOVE_ITEM lib_sources ${PROJECT_SOURCE_DIR}/src/teufel.cpp)
# files not yet finished and not included in the build
list(REMOVE_ITEM lib_sources ${PROJECT_SOURCE_DIR}/src/cavity.cpp)
list(REMOVE_ITEM lib_sources ${PROJECT_SOURCE_DIR}/src/simulation.cpp)
list(REMOVE_ITEM lib_sources ${PROJECT_SOURCE_DIR}/src/wave.cpp)
message(STATUS "Found source files ${lib_sources}")
add_library(${PROJECT_NAME}lib ${lib_sources})

set (TEUFEL_VERSION_MAJOR 0)
set (TEUFEL_VERSION_MINOR 8)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "config.h" @ONLY
  )

add_executable(teufel ${PROJECT_SOURCE_DIR}/src/teufel.cpp)
#                              ${HDF5_CXX_LIBRARY_hdf5}
target_link_libraries(teufel ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIB} ${SDDS_MDBLIB} ${SDDS_MDBCOMMON} ${SDDS_LZMA}
                             ${ZLIB_LIBRARY_RELEASE} )
