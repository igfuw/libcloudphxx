# needed for the OpenMP test to work in C++-only project
# (see http://public.kitware.com/Bug/view.php?id=11910)
cmake_minimum_required(VERSION 2.8.11) # for target_include_directories

project(libcloudph++ CXX C)

# the policies we care about:
# - CMP0025 - make CMake distinguis between Apple and LLVM clang
# - CMP0042 - make CMake use RPATHs on OSX
# - CMP0060 - make CMake always keep absoult RPATHs, even if installing in implicit directory
if(CMAKE_VERSION VERSION_GREATER 2.9)
  cmake_policy(VERSION 3.0)
endif()

list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_MACOSX_RPATH ON) # explicit, since policy CMP0042 didn't work...

# using include() istead of find_package(libcloudph++) to use local CMake code
# and not the system-installed one
#include(${CMAKE_SOURCE_DIR}/libcloudph++-config.cmake)
#if (NOT libcloudphxx_FOUND)
#  message(FATAL_ERROR "local libcloudph++-config.cmake not found!")
#endif()

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

# pointing to local headers
set(CMAKE_CXX_FLAGS "-I${CMAKE_SOURCE_DIR}/include ${CMAKE_CXX_FLAGS}")

# enabling additional warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")

# enabling c++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# config-specific flags
set(CMAKE_CXX_FLAGS_DEBUG " -g -Og -DTHRUST_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO " -Ofast -march=native")
set(CMAKE_CXX_FLAGS_RELEASE " -DNDEBUG -Ofast -march=native -Winline")

#############################################################################################
# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  message(STATUS " OpenMP found")
else()
  message(STATUS "CMAKE cant find OpenMP (known issue for AppleClang on OSX)")
  # maybe some hints are here: https://github.com/src-d/kmcuda/issues/6
  set(OpenMP_CXX_FLAGS " ")
endif()

############################################################################################
# CUDA
find_package(CUDA)
if (NOT CUDA_FOUND)
  message(STATUS "CUDA not found. 

* CUDA support will not be compiled.
* To install CUDA, please try:
*   Debian/Ubuntu: sudo apt-get install nvidia-cuda-toolkit
*   Fedora: TODO
*   Homebrew: TODO
  ")
endif()

set(CUDA_PROPAGATE_HOST_FLAGS OFF) # -std=c++11 will not work, TODO: actually now it would

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
      int DeviceCount;
      cudaGetDeviceCount(&DeviceCount);
      for(int i = 0 ; i < DeviceCount; ++i)
      {
        cudaDeviceProp prop; 
        cudaError_t err = cudaGetDeviceProperties(&prop, i); 
        if (err == cudaSuccess)
          std::cout << prop.major << prop.minor << \" \" ;
        else
        {   
          std::cerr << cudaGetErrorString(err);
          std::cout << \"?\";
        }   
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
  separate_arguments(msg)
  list(REMOVE_ITEM msg "")
  list(REMOVE_DUPLICATES msg)

  # find the lowest compute capability and set it as virtual architecture, see:
  # https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-compilation
  list(SORT msg)
  list(GET msg 0 VIRT)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --gpu-architecture=compute_${VIRT}")

  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --gpu-code=")
  foreach(code ${msg})
    message(STATUS "CUDA capability: ${code}")
    if (code STREQUAL "?")
      set (code "20")
      message(STATUS "CUDA capability check failed, assuming a default of ${code}")
    endif()
    set(arch ${code})
    if (arch EQUAL "21") # see https://devtalk.nvidia.com/default/topic/606195/-solved-nvcc-fatal-unsupported-gpu-architecture-compute_21-/
      set(arch "20")
    endif()
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}sm_${arch},")
  endforeach(code)

  # enabling c++11
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11")
  unset(pfx)
  unset(tmpdir)
  unset(msg)
  unset(status)

  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -I${CMAKE_SOURCE_DIR}/include")
  
  # Release with debug info mode cuda flags
  if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -O3 -use_fast_math
        -Xcompiler=-Ofast,-march=native,-fopenmp")
  # Debug mode cuda flags
  elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g -DTHRUST_DEBUG -lineinfo
        \"-DBOOST_NOINLINE=__attribute__((noinline))\" -Xcompiler=-fopenmp,-g,-Og,-rdynamic")  
  # Release mode cuda flags
  else()
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -DNDEBUG -O3 -use_fast_math
        -Xcompiler=-Ofast,-march=native,-DNDEBUG,-fopenmp")
  endif()

  add_definitions(-DCUDA_FOUND)
endif()

############################################################################################
# Thrust, location of Thrust can be hinted by setting THRUST_INCLUDE_DIR
find_package(Thrust)

if (NOT THRUST_INCLUDE_DIR)
  message(FATAL_ERROR "Thrust library not found. 

* To install Thrust, please try:
*   Debian/Ubuntu: sudo apt-get install libthrust-dev
*   Fedora: TODO
*   Homebrew: TODO
  ")
endif()

# make nvcc include the thrust found here instead of the one shipped with cuda
if (CUDA_FOUND)
  CUDA_INCLUDE_DIRECTORIES("${THRUST_INCLUDE_DIR}")
endif()

############################################################################################
# Boost libraries
find_package(Boost)
if(Boost_FOUND)
#TODO: if boost is not linked in some program, link boost libs to libcloudphxx_lgrngn.so ?
#  set(libcloudphxx_LIBRARIES "${libcloudphxx_LIBRARIES};${Boost_LIBRARIES}")
  set(libcloudphxx_INCLUDE_DIRS "${libcloudphxx_INCLUDE_DIRS};${Boost_INCLUDE_DIRS}")
else()
#TODO: check separately for optional and mandatory components
message(FATAL_ERROR "Boost (or some of its components) not found.

* To insall Boost, please try:
*   Debian/Ubuntu: sudo apt-get install libboost-all-dev
*   Fedora: sudo yum install boost-devel
")
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
  add_custom_target(git_revision.h ALL
   git log -1 "--format=format:#define LIBCLOUDPHXX_GIT_REVISION \"%H\"%n" HEAD > include/libcloudph++/git_revision.h
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
