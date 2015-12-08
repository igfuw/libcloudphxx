cmake_minimum_required(VERSION 2.8.11) # for target_include_directories
project(libcloudph++ CXX C)

# using include() istead of find_package(libcloudph++) to use local CMake code
# and not the system-installed one
include(${CMAKE_SOURCE_DIR}/libcloudph++-config.cmake)
if (NOT libcloudphxx_FOUND)
  message(FATAL_ERROR "local libcloudph++-config.cmake not found!")
endif()

# pointing to local headers
set(CMAKE_CXX_FLAGS "-I${CMAKE_SOURCE_DIR}/include ${CMAKE_CXX_FLAGS}")

# enabling additional warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${libcloudphxx_CXX_FLAGS_RELEASE} -Winline")
  set(CMAKE_CXX_FLAGS_RELEASE "")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "${libcloudphxx_CXX_FLAGS_DEBUG} -DTHRUST_DEBUG")
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
  # trying to figure out compute capability of the current 
  message(STATUS "Trying to obtain CUDA capability of local hardware...")
  set(pfx "compute cabability check")
  execute_process(COMMAND "mktemp" "-d" RESULT_VARIABLE status OUTPUT_VARIABLE tmpdir)
  if (NOT status EQUAL 0) 
    message(FATAL_ERROR "${pfx}: mkdtemp failed")
  endif()
  file(WRITE "${tmpdir}/test.cu" "
    #include <iostream>

    int main() 
    {
      cudaDeviceProp prop; 
      cudaError_t err = cudaGetDeviceProperties(&prop, 0); 
      if (err == cudaSuccess)
        std::cout << prop.major << prop.minor;
      else
        std::cout << \"?\";
    }
  ")
  execute_process(COMMAND "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" "test.cu" WORKING_DIRECTORY ${tmpdir} RESULT_VARIABLE status ERROR_VARIABLE msg)
  if (NOT status EQUAL 0) 
    message(FATAL_ERROR "${pfx}: nvcc failed\n ${msg}")
  endif()
  execute_process(COMMAND "./a.out" WORKING_DIRECTORY ${tmpdir} RESULT_VARIABLE status OUTPUT_VARIABLE msg)
  if (NOT status EQUAL 0)
    message(FATAL_ERROR "${pfx}: test program failed")
  endif()
  message(STATUS "CUDA capability: ${msg}")
  if (msg STREQUAL "?")
    set (msg "20")
    message(STATUS "CUDA capability check failed, assuming a default of ${msg}")
  endif()
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --gpu-code=sm_${msg}")
  if (msg EQUAL "21") # see https://devtalk.nvidia.com/default/topic/606195/-solved-nvcc-fatal-unsupported-gpu-architecture-compute_21-/
    set(msg "20")
  endif()
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --gpu-architecture=compute_${msg}")
  unset(pfx)
  unset(tmpdir)
  unset(msg)
  unset(status)

  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -I${CMAKE_SOURCE_DIR}/include")
  
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -DNDEBUG -O3 -Xcompiler=-Ofast,-march=native,-DNDEBUG") #TODO: fix -use_fast_math, now it causes incorrect Beard velocities
  else()
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g -DTHRUST_DEBUG \"-DBOOST_NOINLINE=__attribute__((noinline))\"")
# do not use -G flag on travis since it causes build to fail; probably due to memory shortage
    if(NOT DEFINED ENV{TRAVIS})
      set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -G")
    endif()
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

install(
  FILES
    libcloudph++-config.cmake
  DESTINATION
    share/libcloudph++
)
