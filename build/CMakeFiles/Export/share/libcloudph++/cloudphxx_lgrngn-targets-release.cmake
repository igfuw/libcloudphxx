#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "clphxx::cloudphxx_lgrngn" for configuration "Release"
set_property(TARGET clphxx::cloudphxx_lgrngn APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(clphxx::cloudphxx_lgrngn PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcloudphxx_lgrngn.so"
  IMPORTED_SONAME_RELEASE "libcloudphxx_lgrngn.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS clphxx::cloudphxx_lgrngn )
list(APPEND _IMPORT_CHECK_FILES_FOR_clphxx::cloudphxx_lgrngn "${_IMPORT_PREFIX}/lib/libcloudphxx_lgrngn.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
