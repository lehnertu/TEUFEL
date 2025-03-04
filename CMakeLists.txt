cmake_minimum_required (VERSION 2.6)
project (TEUFEL)

set (TEUFEL_VERSION_MAJOR 0)
set (TEUFEL_VERSION_MINOR 11)
set (TEUFEL_VERSION_PATCH 1)

# where to find modules
# we add teufel/lib to the search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/lib")
message(STATUS "cmake search path: ${CMAKE_MODULE_PATH}")

set(CMAKE_VERBOSE_MAKEFILE OFF)

# Set default build type.
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "CMAKE_BUILD_TYPE not given; setting to 'Release'.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# HDF5 library
find_package(HDF5 COMPONENTS CXX REQUIRED)
if (HDF5_FOUND)
    message(STATUS "HDF5 includes ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5 libraries ${HDF5_LIBRARIES}")
endif (HDF5_FOUND)
find_package(ZLIB REQUIRED)

# XML library
find_package(PUGIXML REQUIRED)

# muParser library
find_package(MUPARSER REQUIRED)

# our local SDDS library
set(SDDS_DIR ${PROJECT_SOURCE_DIR}/lib/SDDS)
find_library(SDDS_LIB NAMES SDDS1 HINTS ${SDDS_DIR})
find_library(SDDS_MDBLIB NAMES mdblib HINTS ${SDDS_DIR})
find_library(SDDS_MDBCOMMON NAMES mdbcommon HINTS ${SDDS_DIR})
find_library(SDDS_LZMA NAMES lzma HINTS ${SDDS_DIR})
set (SDDS_LIBRARIES ${SDDS_LIB} ${SDDS_MDBLIB} ${SDDS_MDBCOMMON} ${SDDS_LZMA})
find_path(SDDS_INCLUDE_DIR NAMES SDDS.h HINTS ${SDDS_DIR})
if(SDDS_INCLUDE_DIR)
    message(STATUS "SDDS includes found in ${SDDS_INCLUDE_DIR}")
endif()

# EIGEN library
find_package (Eigen3 3.3 REQUIRED)
find_package (Eigen3)
if (EIGEN3_FOUND)
    message(STATUS "Eigen3 includes ${EIGEN3_INCLUDE_DIR}")
else (EIGEN3_FOUND)
    message(STATUS "Eigen3 not found")
endif (EIGEN3_FOUND)

# Open-MP
include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(OPENMP_FOUND)
  message("ERROR: OpenMP could not be found.")
endif(OPENMP_FOUND)

# MPI library
find_package(MPI COMPONENTS CXX REQUIRED)
if (MPI_CXX_FOUND)
    message(STATUS "MPI includes ${MPI_CXX_INCLUDE_PATH}")
    message(STATUS "MPI libraries ${MPI_CXX_LIBRARIES}")
    message(STATUS "MPI compile flags ${MPI_CXX_COMPILE_FLAGS}")
    message(STATUS "MPI link flags ${MPI_CXX_LINK_FLAGS}")
endif ()

# where to find header files
include_directories( ${PROJECT_BINARY_DIR}
                     ${PROJECT_SOURCE_DIR}/src
                     ${EIGEN3_INCLUDE_DIR}
                     ${HDF5_INCLUDE_DIRS}
                     ${PUGIXML_INCLUDE_DIR}
                     ${MUPARSER_INCLUDE_DIR}
                     ${SDDS_INCLUDE_DIR}
                     ${MPI_CXX_INCLUDE_PATH})

# we put all objects into a library that can be linked to the various executables
aux_source_directory(${PROJECT_SOURCE_DIR}/src lib_sources )
# the main program is not part of the library
list(REMOVE_ITEM lib_sources ${PROJECT_SOURCE_DIR}/src/teufel_hybrid.cpp)
# stand-alone codes are not part of the library
list(REMOVE_ITEM lib_sources ${PROJECT_SOURCE_DIR}/src/propagate.cpp)
list(REMOVE_ITEM lib_sources ${PROJECT_SOURCE_DIR}/src/reflection.cpp)
# files not yet finished and not included in the build
list(REMOVE_ITEM lib_sources ${PROJECT_SOURCE_DIR}/src/cavity.cpp)
message(STATUS "Found source files ${lib_sources}")
add_library(${PROJECT_NAME}lib ${lib_sources})

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "config.h" @ONLY
  )

# compiler flags
add_definitions(-std=c++11 -Wall -m64 -pthread)

add_executable(teufel ${PROJECT_SOURCE_DIR}/src/teufel_hybrid.cpp)
set_target_properties(teufel PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set_target_properties(teufel PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
target_link_libraries(teufel
                             ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             ${MPI_CXX_LIBRARIES}
                             ${PUGIXML_LIBRARIES}
                             ${MUPARSER_LIBRARY}
                             ${ZLIB_LIBRARY_RELEASE} )

add_executable(propagate ${PROJECT_SOURCE_DIR}/src/propagate.cpp)
set_target_properties(propagate PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set_target_properties(propagate PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
target_link_libraries(propagate ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             Eigen3::Eigen
                             ${PUGIXML_LIBRARIES}
                             ${MUPARSER_LIBRARY}
                             ${ZLIB_LIBRARY_RELEASE} )

add_executable(reflection ${PROJECT_SOURCE_DIR}/src/reflection.cpp)
set_target_properties(reflection PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set_target_properties(reflection PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
target_link_libraries(reflection ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             Eigen3::Eigen
                             ${PUGIXML_LIBRARIES}
                             ${MUPARSER_LIBRARY}
                             ${ZLIB_LIBRARY_RELEASE} )

# for building the documentation
option(BUILD_DOCUMENTATION "Generate makefiles for documentation" ON)
option(BUILD_DOCUMENTATION_DEFAULT "Build documentation by default" OFF)
if(BUILD_DOCUMENTATION)
    find_package(Doxygen)
        if(NOT DOXYGEN_FOUND)
            message(FATAL_ERROR "Doxygen is not installed or not properly configured")
        endif(NOT DOXYGEN_FOUND)
    find_package(LATEX COMPONENTS PDFLATEX)
        if(NOT LATEX_FOUND)
            message(FATAL_ERROR "LaTeX is not installed or not properly configured")
        endif(NOT LATEX_FOUND)
    set(DOC_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doc)
    message(DOC_SRC_DIR = ${DOC_SRC_DIR})
    set(DOC_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/doc)
    message(DOC_BUILD_DIR = ${DOC_BUILD_DIR})
    file(MAKE_DIRECTORY ${DOC_BUILD_DIR})
    if(BUILD_DOCUMENTATION_DEFAULT)
        add_custom_target(doc ALL
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/setup.dox
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen")
        add_custom_target(man ALL
            # we first compile with LaTeX to generate all indices
            COMMAND ${LATEX_COMPILER} ${DOC_SRC_DIR}/UserManual.tex 
            # then we use PDFLATEX to generate the final output
            COMMAND ${PDFLATEX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/doc/UserManual.tex
            COMMAND cp -f UserManual.pdf ${DOC_SRC_DIR}/
            WORKING_DIRECTORY ${DOC_BUILD_DIR}
            COMMENT "Generating user manual with LaTeX"
            )
    else()
        add_custom_target(doc
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/setup.dox
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen")
        add_custom_target(man
            # we first compile with LaTeX to generate all indices
            COMMAND ${LATEX_COMPILER} ${DOC_SRC_DIR}/UserManual.tex 
            # then we use PDFLATEX to generate the final output
            COMMAND ${PDFLATEX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/doc/UserManual.tex
            COMMAND cp -f UserManual.pdf ${DOC_SRC_DIR}/
            WORKING_DIRECTORY ${DOC_BUILD_DIR}
            COMMENT "Generating user manual with LaTeX"
            )
    endif()
endif()

# for building test binaries
option(BUILD_TESTS "Build binaries for the test cases" ON)
if(BUILD_TESTS)
    set(PROJECT_TESTS_DIR ${PROJECT_SOURCE_DIR}/tests)
    add_executable(teufel.integrate_field ${PROJECT_TESTS_DIR}/teufel.integrate_field.cpp)
    target_link_libraries(teufel.integrate_field ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.bunch ${PROJECT_TESTS_DIR}/teufel.bunch.cpp)
    target_link_libraries(teufel.bunch ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.magnet ${PROJECT_TESTS_DIR}/teufel.magnet.cpp)
    target_link_libraries(teufel.magnet ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.undulator ${PROJECT_TESTS_DIR}/teufel.undulator.cpp)
    target_link_libraries(teufel.undulator ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             ${PUGIXML_LIBRARIES}
                             ${MUPARSER_LIBRARY}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.loop ${PROJECT_TESTS_DIR}/teufel.loop.cpp)
    target_link_libraries(teufel.loop ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.EcrossB ${PROJECT_TESTS_DIR}/teufel.EcrossB.cpp)
    target_link_libraries(teufel.EcrossB ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.RadPush ${PROJECT_TESTS_DIR}/teufel.RadPush.cpp)
    target_link_libraries(teufel.RadPush ${PROJECT_NAME}lib
                             ${PUGIXML_LIBRARIES}
                             ${MUPARSER_LIBRARY}
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.electrostatic ${PROJECT_TESTS_DIR}/teufel.electrostatic.cpp)
    target_link_libraries(teufel.electrostatic ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             ${ZLIB_LIBRARY_RELEASE} )
    add_executable(teufel.diffraction ${PROJECT_TESTS_DIR}/teufel.DiffractionScreen.cpp)
    target_link_libraries(teufel.diffraction ${PROJECT_NAME}lib
                             ${HDF5_LIBRARIES}
                             ${SDDS_LIBRARIES}
                             ${ZLIB_LIBRARY_RELEASE} )
endif()
