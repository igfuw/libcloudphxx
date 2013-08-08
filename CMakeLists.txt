cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(libcloudph++ CXX C)

# requireing C++11 and pointing to local headers
set(CMAKE_CXX_FLAGS "-std=c++11 -I${CMAKE_SOURCE_DIR}/include ${CMAKE_CXX_FLAGS}")

# TODO: libmpdataxx-config?
# -DNDEBUG?
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_RELEASE "-Wfatal-errors -Ofast -march=native")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "-Wfatal-errors -O3 -march=native -ffast-math") 
endif()


add_subdirectory(include) 
add_subdirectory(src) 

#TODO: search for Blitz
#TODO: search for libmpdata++


# TODO: automate "make doc"

# tests
enable_testing()
add_subdirectory(tests)
