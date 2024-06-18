#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "simfor::simfor" for configuration ""
set_property(TARGET simfor::simfor APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(simfor::simfor PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libsimfor.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS simfor::simfor )
list(APPEND _IMPORT_CHECK_FILES_FOR_simfor::simfor "${_IMPORT_PREFIX}/lib/libsimfor.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
