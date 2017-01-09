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

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS_DEBUG "${libcloudphxx_CXX_FLAGS_DEBUG} -DTHRUST_DEBUG")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${libcloudphxx_CXX_FLAGS_RELWITHDEBINFO}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${libcloudphxx_CXX_FLAGS_RELEASE} -Winline")
  set(CMAKE_CXX_FLAGS_RELEASE "")
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
#find_package(CUDA QUIET)
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
  execute_process(COMMAND "mktemp" "-d" "/tmp/tmp.XXX" RESULT_VARIABLE status OUTPUT_VARIABLE tmpdir)
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
      {
        std::cerr << cudaGetErrorString(err);
        std::cout << \"?\";
      }
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

  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -I${CMAKE_SOURCE_DIR}/include -std=c++11")
  
  # Release with debug info mode cuda flags
  if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -O3 -use_fast_math
        -Xcompiler=-Ofast,-march=native,-fopenmp")
  # Debug mode cuda flags
  elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g -DTHRUST_DEBUG -lineinfo
        \"-DBOOST_NOINLINE=__attribute__((noinline))\" -Xcompiler=-fopenmp,-g,-rdynamic")  
  # Release mode cuda flags
  else()
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -DNDEBUG -O3 -use_fast_math
        -Xcompiler=-Ofast,-march=native,-DNDEBUG,-fopenmp")
  endif()

  add_definitions(-DCUDA_FOUND)
endif()

############################################################################################
# BOOST ODEINT VERSION TEST
message(STATUS "Testing if Boost ODEINT version >= 1.58")
set(pfx "boost odeint check")
execute_process(COMMAND "mktemp" "-d" "/tmp/tmp.XXX" RESULT_VARIABLE status OUTPUT_VARIABLE tmpdir)
if (NOT status EQUAL 0) 
  message(FATAL_ERROR "${pfx}: mkdtemp failed")
endif()
file(WRITE "${tmpdir}/test.cpp" "
  #define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CPP
  
  #include <thrust/system/cpp/vector.h>

  #include <boost/numeric/odeint.hpp>
  #include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
  #include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
  #include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>
  
  struct rhs 
  {
    void operator()(
      const thrust::cpp::vector<float> &psi,
      thrust::cpp::vector<float> &dot_psi,
      const float /* t */
    )  
    {
    }  
  };
  
  int main() 
  { 
    boost::numeric::odeint::euler<
      thrust::cpp::vector<float>, // state_type
      float,                      // value_type
      thrust::cpp::vector<float>, // deriv_type
      float,                      // time_type
      boost::numeric::odeint::thrust_algebra,
      boost::numeric::odeint::thrust_operations 
    > chem_stepper;
  
    thrust::cpp::vector<float> v(11);
    chem_stepper.do_step(rhs(), v, 0, 1);
  }
")
execute_process(COMMAND "${CMAKE_CXX_COMPILER}" "test.cpp" "-I${Boost_INCLUDE_DIR}" WORKING_DIRECTORY ${tmpdir} RESULT_VARIABLE status ERROR_VARIABLE msg)
if (NOT status EQUAL 0) 
  message(FATAL_ERROR "${pfx}: c++ compiler failed\n ${msg}")
endif()
execute_process(COMMAND "./a.out" WORKING_DIRECTORY ${tmpdir} RESULT_VARIABLE status OUTPUT_VARIABLE msg)
if (NOT status EQUAL 0)
  message(FATAL_ERROR "${pfx}: test program failed, install Boost odeint version >= 1.58")
endif()
unset(pfx)
unset(tmpdir)
unset(msg)
unset(status)

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
