# Find the muparser equation parsing library.
# https://github.com/beltoforion/muparser.git
#
# Sets the usual variables expected for find_package scripts:
#
# MUPARSER_INCLUDE_DIR - header location
# MUPARSER_LIBRARIES - library to link against
# MUPARSER_FOUND - true if muparser was found.

# find_package(PkgConfig)
# pkg_check_modules(MUPARSER REQUIRED muparser)
# find_path(MUPARSER_INCLUDE_DIR muParser.h
#           HINTS ${PMUPARSER_INCLUDEDIR} ${MUPARSER_INCLUDE_DIRS})
# find_library(MUPARSER_LIBRARY NAMES muparser
#           HINTS ${MUPARSER_LIBDIR} ${MUPARSER_LIBRARY_DIRS})
# message(STATUS "muParser include : ${MUPARSER_INCLUDE_DIR}")
# message(STATUS "muParser library : ${MUPARSER_LIBRARY}")

message(STATUS "this is FindMUPARSER.cmake")

unset (MUPARSER_LIBRARY CACHE)
unset (MUPARSER_INCLUDE_DIR CACHE)

find_path (MUPARSER_INCLUDE_DIR
           NAMES muParser.h
           PATHS /usr/include
           /usr/local/include
           ${CMAKE_SOURCE_DIR}/lib/muparser/include)
message(STATUS "MUPARSER_INCLUDE_DIR = ${MUPARSER_INCLUDE_DIR}")
           
find_library (MUPARSER_LIBRARY
              NAMES libmuparser.a libmuparser.so
              PATHS /usr/lib
              /usr/lib/x86_64-linux-gnu/
              /usr/local/lib
              ${CMAKE_SOURCE_DIR}/lib/muparser/build)
message(STATUS "MUPARSER_LIBRARY = ${MUPARSER_LIBRARY}")

# Support the REQUIRED and QUIET arguments, and set MUPARSER_FOUND if found.
include (FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS (MUPARSER DEFAULT_MSG MUPARSER_LIBRARY
                                   MUPARSER_INCLUDE_DIR)

if (MUPARSER_FOUND)
    set (MUPARSER_LIBRARIES ${MUPARSER_LIBRARY})
    if (NOT MUPARSER_FIND_QUIETLY)
        message (STATUS "MUPARSER include = ${MUPARSER_INCLUDE_DIR}")
        message (STATUS "MUPARSER library = ${MUPARSER_LIBRARY}")
    endif ()
else ()
    message (STATUS "No MUPARSER found")
endif()

mark_as_advanced (MUPARSER_LIBRARY MUPARSER_INCLUDE_DIR)

