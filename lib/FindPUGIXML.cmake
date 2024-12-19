# Find the pugixml XML parsing library.
#
# Sets the usual variables expected for find_package scripts:
#
# PUGIXML_INCLUDE_DIR - header location
# PUGIXML_LIBRARIES - library to link against
# PUGIXML_FOUND - true if pugixml was found.

message(STATUS "this is FindPUGIXML.cmake")

unset (PUGIXML_LIBRARY CACHE)
unset (PUGIXML_INCLUDE_DIR CACHE)

# If PUGIXML_ROOT is not yet known but was defined in the environment, use that.
if (NOT PUGIXML_ROOT AND NOT $ENV{PUGIXML_ROOT} STREQUAL "")
    SET(PUGIXML_ROOT $ENV{PUGIXML_ROOT})
endif()

find_path (PUGIXML_INCLUDE_DIR
           NAMES pugixml.hpp
           PATHS ${PUGIXML_ROOT}
           PATHS ${PUGIXML_ROOT}/include
           /usr/include
           /usr/local/include
           ${CMAKE_SOURCE_DIR}/lib/pugixml)
           
find_library (PUGIXML_LIBRARY
              NAMES pugixml
              PATHS ${PUGIXML_ROOT}
              PATHS ${PUGIXML_ROOT}/lib
              /usr/lib
              /usr/lib/x86_64-linux-gnu/
              /usr/local/lib
              ${CMAKE_SOURCE_DIR}/lib/pugixml)

# Support the REQUIRED and QUIET arguments, and set PUGIXML_FOUND if found.
include (FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS (PUGIXML DEFAULT_MSG PUGIXML_LIBRARY
                                   PUGIXML_INCLUDE_DIR)

if (PUGIXML_FOUND)
    set (PUGIXML_LIBRARIES ${PUGIXML_LIBRARY})
    if (NOT PUGIXML_FIND_QUIETLY)
        message (STATUS "PugiXML include = ${PUGIXML_INCLUDE_DIR}")
        message (STATUS "PugiXML library = ${PUGIXML_LIBRARY}")
    endif ()
else ()
    message (STATUS "No PugiXML found")
endif()

mark_as_advanced (PUGIXML_LIBRARY PUGIXML_INCLUDE_DIR)

