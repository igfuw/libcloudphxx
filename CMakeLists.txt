cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(libcloudph++ CXX C)

# requireing C++11 
set(CMAKE_CXX_FLAGS "-std=c++11")

add_subdirectory(include) 
add_subdirectory(src) 

#TODO: search for Blitz
#TODO: search for libmpdata++


# TODO: automate "make doc"

# tests
enable_testing()
add_subdirectory(tests)
