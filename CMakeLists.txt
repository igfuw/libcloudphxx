cmake_minimum_required(VERSION 2.8.11) # for target_include_directories
project(libcloudph++ CXX C)

# requireing C++11 and pointing to local headers
set(CMAKE_CXX_FLAGS "-std=c++11 -Wfatal-errors -I${CMAKE_SOURCE_DIR}/include ${CMAKE_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -march=native")
endif()

find_package(CUDA QUIET)
set(CUDA_PROPAGATE_HOST_FLAGS OFF) # -std=c++11 will not work
# TODO: propagate optimisation flags to nvcc

# TODO: detecd CUDA arch (native functionality would be best?)
if (CUDA_FOUND)
  set(CUDA_NVCC_FLAGS "-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -I${CMAKE_SOURCE_DIR}/include ${CUDA_NVCC_FLAGS}")
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -DNDEBUG")
  endif()
  add_definitions(-DCUDA_FOUND)
endif()

add_subdirectory(src) 

enable_testing()
add_subdirectory(tests)

add_subdirectory(include) 
