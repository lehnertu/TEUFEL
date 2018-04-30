cmake_minimum_required (VERSION 2.6)
project (TEUFEL)

set(CMAKE_VERBOSE_MAKEFILE ON)

# Set default build type.
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "CMAKE_BUILD_TYPE not given; setting to 'Release'.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# for building the documentation
option(BUILD_DOCUMENTATION "Generate doxygen documentation" OFF)
find_package(Doxygen)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is not installed or not properly configured")
    endif(NOT DOXYGEN_FOUND)
if(BUILD_DOCUMENTATION)
    add_custom_target(docs ALL
                      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/setup.dox
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                      COMMENT "Generating API documentation with Doxygen")
else()
    add_custom_target(docs
                      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/setup.dox
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                      COMMENT "Generating API documentation with Doxygen")
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
target_link_libraries(teufel ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIB} ${SDDS_MDBLIB} ${SDDS_MDBCOMMON} ${SDDS_LZMA}
                             ${ZLIB_LIBRARY_RELEASE} )

option(BUILD_TESTS "Build binaries for the test cases" ON)
if(BUILD_TESTS)
    set(PROJECT_TESTS_DIR ${PROJECT_SOURCE_DIR}/tests)
    add_executable(teufel.integrate_field ${PROJECT_TESTS_DIR}/teufel.integrate_field.cpp)
    target_link_libraries(teufel.integrate_field ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIB} ${SDDS_MDBLIB} ${SDDS_MDBCOMMON} ${SDDS_LZMA}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.bunch ${PROJECT_TESTS_DIR}/teufel.bunch.cpp)
    target_link_libraries(teufel.bunch ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIB} ${SDDS_MDBLIB} ${SDDS_MDBCOMMON} ${SDDS_LZMA}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.magnet ${PROJECT_TESTS_DIR}/teufel.magnet.cpp)
    target_link_libraries(teufel.magnet ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIB} ${SDDS_MDBLIB} ${SDDS_MDBCOMMON} ${SDDS_LZMA}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.undulator ${PROJECT_TESTS_DIR}/teufel.undulator.cpp)
    target_link_libraries(teufel.undulator ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIB} ${SDDS_MDBLIB} ${SDDS_MDBCOMMON} ${SDDS_LZMA}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.loop ${PROJECT_TESTS_DIR}/teufel.loop.cpp)
    target_link_libraries(teufel.loop ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIB} ${SDDS_MDBLIB} ${SDDS_MDBCOMMON} ${SDDS_LZMA}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.EcrossB ${PROJECT_TESTS_DIR}/teufel.EcrossB.cpp)
    target_link_libraries(teufel.EcrossB ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIB} ${SDDS_MDBLIB} ${SDDS_MDBCOMMON} ${SDDS_LZMA}
                             ${ZLIB_LIBRARY_RELEASE} )
endif()

option(BUILD_EXAMPLES "Build binaries for the examples" ON)
if(BUILD_EXAMPLES)
else()
endif()
