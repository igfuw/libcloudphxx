#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "clphxx::cloudphxx_lgrngn" for configuration "Debug"
set_property(TARGET clphxx::cloudphxx_lgrngn APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(clphxx::cloudphxx_lgrngn PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libcloudphxx_lgrngn_dbg.so"
  IMPORTED_SONAME_DEBUG "libcloudphxx_lgrngn_dbg.so"
  )

list(APPEND _cmake_import_check_targets clphxx::cloudphxx_lgrngn )
list(APPEND _cmake_import_check_files_for_clphxx::cloudphxx_lgrngn "${_IMPORT_PREFIX}/lib/libcloudphxx_lgrngn_dbg.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
