cmake_minimum_required(VERSION 3.0)
project(icicle CXX C)
enable_testing()

find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
    message(FATAL_ERROR "libclouphxx++ not found")
endif()

add_subdirectory(tests)

