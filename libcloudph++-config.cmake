
############################################################################################
# the following variables will be set:
set(libcloudphxx_FOUND False)
set(libcloudphxx_INCLUDE_DIRS "")
set(libcloudphxx_LIBRARIES "")

############################################################################################
# libcloudphxx libs and headers 
# also work for non-default install location (i.e. for make DESTDIR=<dir> install)
set(libcloudphxx_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../include/")
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CONFIG_SUFFIX "")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(CONFIG_SUFFIX "_relwithdbg")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CONFIG_SUFFIX "_dbg")
endif()

set(libcloudphxx_LIBRARIES "${CMAKE_CURRENT_LIST_DIR}/../../lib/libcloudphxx_lgrngn${CONFIG_SUFFIX}.so")
if(NOT EXISTS libcloudphxx_LIBRARIES)
  message(FATAL_ERROR "The libcloudph++ library for selected config not found at ${libcloudphxx_LIBRARIES}") 
endif() 


############################################################################################
# Boost libraries
find_package(Boost)
if(Boost_FOUND)
#TODO: if boost is not linked in some program, link boost libs to libcloudphxx_lgrngn.so ?
#  set(libcloudphxx_LIBRARIES "${libcloudphxx_LIBRARIES};${Boost_LIBRARIES}")
  set(libcloudphxx_INCLUDE_DIRS "${libcloudphxx_INCLUDE_DIRS};${Boost_INCLUDE_DIRS}")
else()
#TODO: check separately for optional and mandatory components
message(FATAL_ERROR "Boost (or some of its components) not found.

* To insall Boost, please try:
*   Debian/Ubuntu: sudo apt-get install libboost-all-dev
*   Fedora: sudo yum install boost-devel
")
endif()

############################################################################################
#list(REMOVE_DUPLICATES libcloudpxxx_INCLUDE_DIRS)
list(REMOVE_ITEM libcloudphxx_INCLUDE_DIRS "")
set(libcloudphxx_FOUND TRUE)
