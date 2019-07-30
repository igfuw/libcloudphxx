cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(icicle CXX C)
enable_testing()

# search for a local version of libcloudphxx
# TODO: do it using the libcloudph++ target? what about proper lib for selected condig type?
#unset(libcloudphxx_LIBRARIES CACHE)
#set(libcloud_postfix "")
#if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
#  set(libcloud_postfix "_relwithdbg")
#elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
#  set(libcloud_postfix "_dbg")
#endif()
#find_library(libcloudphxx_LIBRARIES cloudphxx_lgrngn${libcloud_postfix} ${CMAKE_SOURCE_DIR}/../../build NO_DEFAULT_PATH)
#if (libcloudphxx_LIBRARIES)
#  message(STATUS "icicle will use libclouphxx++ from current repository")
#  include_directories(PUBLIC ${CMAKE_SOURCE_DIR}/../../include)
#else()
#  message(STATUS "libcloudphxx++ not found in current repository, icicle will use an installed version")
#endif()

# if no local version found, look for an installed one
#find_library(libcloudphxx_LIBRARIES cloudphxx_lgrngn)

find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

find_package(libmpdata++)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  # make the Release flags the default
  set(CMAKE_CXX_FLAGS ${libmpdataxx_CXX_FLAGS_RELEASE})
  # ignore CMake defaults
  set(CMAKE_CXX_FLAGS_RELEASE "")
else()
  set(CMAKE_CXX_FLAGS_DEBUG ${libmpdataxx_CXX_FLAGS_DEBUG})
endif()

add_subdirectory(src)
add_subdirectory(tests)
