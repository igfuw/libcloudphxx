cmake_minimum_required(VERSION 2.8.11) # for target_include_directories
project(libcloudph++ CXX C)

# using include() istead of find_package(libcloudph++) to use local CMake code
# and not the system-installed one
include(${CMAKE_SOURCE_DIR}/libcloudph++-config.cmake)
if (NOT libcloudphxx_FOUND)
  message(FATAL_ERROR "local libcloudph++-config.cmake not found!")
endif()

# requireing C++11 and pointing to local headers
set(CMAKE_CXX_FLAGS "-I${CMAKE_SOURCE_DIR}/include ${CMAKE_CXX_FLAGS}")

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${libcloudphxx_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELEASE "")
else()
  set(CMAKE_CXX_FLAGS_DEBUG ${libcloudphxx_CXX_FLAGS_DEBUG})
endif()

############################################################################################
# Thrust
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
    #define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CPP
    #include <thrust/version.h>
    int main() {}
  " THRUST_FOUND)
if (NOT THRUST_FOUND)
  message(FATAL_ERROR "Thrust library not found. 

* To install Thrust, please try:
*   Debian/Ubuntu: sudo apt-get install libthrust-dev
*   Fedora: TODO
*   Homebrew: TODO
  ")
endif()

############################################################################################
# CUDA
find_package(CUDA QUIET)
if (NOT CUDA_FOUND)
  message(STATUS "CUDA not found. 

* CUDA support will not be compiled.
* To install CUDA, please try:
*   Debian/Ubuntu: sudo apt-get install nvidia-cuda-toolkit
*   Fedora: TODO
*   Homebrew: TODO
  ")
endif()

set(CUDA_PROPAGATE_HOST_FLAGS OFF) # -std=c++11 will not work

if (CUDA_FOUND)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_20,code=sm_20")
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_30,code=sm_30")
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -I${CMAKE_SOURCE_DIR}/include")
  
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -DNDEBUG -O3 -use_fast_math -Xcompiler=\"-Ofast -march=native -DNDEBUG\"")
  else()
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g \"-DBOOST_NOINLINE=__attribute__((noinline))\"")
  endif()

  add_definitions(-DCUDA_FOUND)
endif()

# generate a header file with git revision id
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
  add_custom_target(git_revision.h
   git log -1 "--format=format:#define GIT_REVISION \"%H\"%n" HEAD > include/libcloudph++/git_revision.h
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} VERBATIM
  )
endif()

add_subdirectory(src) 

enable_testing()
add_subdirectory(tests)

add_subdirectory(include) 
add_subdirectory(bindings) 
