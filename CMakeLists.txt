cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(libcloudph++ CXX C)

# requireing C++11 and pointing to local headers
set(CMAKE_CXX_FLAGS "-std=c++11 -Wfatal-errors -I${CMAKE_SOURCE_DIR}/include ${CMAKE_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -march=native")
endif()

add_subdirectory(src) 

enable_testing()
add_subdirectory(tests)

add_subdirectory(include) 

#TODO: search for libmpdata++
