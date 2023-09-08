cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(icicle CXX C)
enable_testing()

find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
    message(FATAL_ERROR "libclouphxx++ not found")
endif()

add_subdirectory(src)
add_subdirectory(tests)

