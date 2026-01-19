# needed for the OpenMP test to work in C++-only project
# (see http://public.kitware.com/Bug/view.php?id=11910)
cmake_minimum_required(VERSION 3.13 FATAL_ERROR) # 3.13 for CMP0079

project(libcloudph++ LANGUAGES CXX)

# the policies we care about:
# - CMP0025 - make CMake distinguis between Apple and LLVM clang
# - CMP0042 - make CMake use RPATHs on OSX
# - CMP0060 - make CMake always keep absoult RPATHs, even if installing in implicit directory
# - CMP0079 - make CMake allow target_link_libraries on targets created in other directories
#if(CMAKE_VERSION VERSION_GREATER 2.9)
#  cmake_policy(VERSION 3.0)
#endif()

cmake_policy(SET CMP0079 NEW)

# put out custom find scripts first
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)

#set(CMAKE_MACOSX_RPATH ON) # explicit, since policy CMP0042 didn't work...
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

# clearing default config-specific flags
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFOPORTABLE "") # without march=native

############################################################################################
# see if CUDA is available

option(LIBCLOUDPHXX_DISABLE_CUDA "LIBCLOUDPHXX_DISABLE_CUDA" False)

if (NOT LIBCLOUDPHXX_DISABLE_CUDA)
  include(CheckLanguage)
  check_language(CUDA)
endif()

if (NOT CMAKE_CUDA_COMPILER)
  message(STATUS "CUDA not found. 

* CUDA support will not be compiled.
* To install CUDA, please try:
*   Debian/Ubuntu: sudo apt-get install nvidia-cuda-toolkit
*   Fedora: TODO
*   Homebrew: TODO
  ")
else()
  enable_language(CUDA)
  # Check CUDA version
  # if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 12.2)
  #     message(FATAL_ERROR "CUDA >= 12.2 is required, but found ${CMAKE_CUDA_COMPILER_VERSION}") # needed for cuda::std::identity
  # endif()
  set(CMAKE_CUDA_FLAGS_DEBUG "")
  set(CMAKE_CUDA_FLAGS_RELEASE "")
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "")
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFOPORTABLE "")
endif()

#############################################################################################
# see if OpenMP is available

find_package(OpenMP)
if(OPENMP_CXX_FOUND)
  message(STATUS " OpenMP found")
else()
  message(STATUS "CMAKE cant find OpenMP (known issue for AppleClang on OSX)")
  # maybe some hints are here: https://github.com/src-d/kmcuda/issues/6
  set(OpenMP_CXX_FLAGS " ")
endif()

############################################################################################
# allowing runtime choice between CUDA, CPP and OpenMP backends

set(files "")

set(files "${files};src/lib.cpp")

set(files "${files};src/lib_cpp.cpp")

if (OPENMP_CXX_FOUND)
  set(files "${files};src/lib_omp.cpp")
endif()

if (CMAKE_CUDA_COMPILER)
  set(files "${files};src/lib_cuda.cu")
endif()

set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER}) # mpic++ will not work

############################################################################################
# create the lagrangian microphysics library target
add_library(cloudphxx_lgrngn SHARED ${files}) 

############################################################################################

# pointing to library headers
target_include_directories(cloudphxx_lgrngn 
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

# enabling additional warnings
#target_compile_options(cloudphxx_lgrngn PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -Wextra>)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")

# enabling c++14
target_compile_features(cloudphxx_lgrngn PUBLIC cxx_std_14)
# config-specific flags
target_compile_options(cloudphxx_lgrngn PRIVATE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:                  -g -Og -DTHRUST_DEBUG>)
target_compile_options(cloudphxx_lgrngn PRIVATE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:                -Ofast -march=native -DNDEBUG >)
target_compile_options(cloudphxx_lgrngn PRIVATE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFO>>:         -Ofast -march=native>)
target_compile_options(cloudphxx_lgrngn PRIVATE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFOPORTABLE>>: -Ofast>)

# add dependencies
if (CMAKE_CUDA_COMPILER)
  target_link_libraries(cloudphxx_lgrngn PRIVATE curand)
endif()

if(OPENMP_CXX_FOUND)
  # target_link_libraries would be preferable, but it adds -fopenmp to nvcc and that breaks the build
  #target_link_libraries(cloudphxx_lgrngn PRIVATE OpenMP::OpenMP_CXX>)

  target_compile_options(cloudphxx_lgrngn PUBLIC $<$<COMPILE_LANGUAGE:CXX>: ${OpenMP_CXX_FLAGS}>)
  target_link_libraries(cloudphxx_lgrngn PUBLIC ${OpenMP_CXX_LIBRARIES})
endif()

############################################################################################
# MPI - detecting if the C++ compiler is actually an MPI wrapper

set(msg "Detecting if the compiler is an MPI wrapper...")
message(STATUS "${msg}")
execute_process(COMMAND ${CMAKE_CXX_COMPILER} "-show" RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_QUIET)
if (status EQUAL 0 AND output MATCHES "mpi")
  set(USE_MPI TRUE)
  target_compile_definitions(cloudphxx_lgrngn PRIVATE USE_MPI)
else()
  set(USE_MPI FALSE)
endif()
message(STATUS "${msg} - ${USE_MPI}")
unset(msg)
unset(status)
unset(output)

############################################################################################
# figure out CUDA nvcc flags
if (CMAKE_CUDA_COMPILER)
  # trying to figure out compute capability of the current 
  message(STATUS "Trying to obtain CUDA capability of local hardware...")
  set(pfx "compute cabability check")
  execute_process(COMMAND "mktemp" "-d" "/tmp/tmp.XXX" RESULT_VARIABLE status OUTPUT_VARIABLE tmpdir)
  string(STRIP "${tmpdir}" tmpdir)  # Removes leading/trailing whitespace and newlines
  if (NOT status EQUAL 0) 
    message(FATAL_ERROR "${pfx}: mkdtemp failed")
  endif()
  file(WRITE "${tmpdir}/test.cu" "
    #include <iostream>
    #include <fstream>
    
    int main() 
    {
      int DeviceCount;
      std::ofstream fdev(\"${tmpdir}/device_capability.txt\");
      cudaError_t error_id = cudaGetDeviceCount(&DeviceCount);

      if (error_id != cudaSuccess) { // error_id==100 is 'no CUDA-capable device is detected'
        std::cerr << \"cudaGetDeviceCount returned \" <<  static_cast<int>(error_id) << \" -> \" << cudaGetErrorString(error_id) << std::endl;
        return 0;
  //      exit(EXIT_FAILURE);
      }

      int driverVersion, runtimeVersion;
      for(int i = 0 ; i < DeviceCount; ++i)
      {
        cudaSetDevice(i);
        cudaDeviceProp prop; 
        cudaError_t error_id = cudaGetDeviceProperties(&prop, i); 
        if (error_id == cudaSuccess)
        {
          std::cout << \"    Device \" << i << \": \" << prop.name << std::endl;
          cudaDriverGetVersion(&driverVersion);
          cudaRuntimeGetVersion(&runtimeVersion);
          std::cout << \"      CUDA Driver Version / Runtime Version  \" 
            << driverVersion / 1000 << \".\" << (driverVersion % 100) / 10 << \" / \" 
            <<  runtimeVersion / 1000 << \".\" << (runtimeVersion % 100) / 10 << std::endl;
          std::cout << \"      CUDA Capability Major/Minor version number: \" << prop.major << \".\" << prop.minor << std::endl;
          fdev << prop.major << prop.minor << \" \" ;
        }
        else
        {   
          std::cerr << \"cudaGetDeviceProperties returned \" <<  static_cast<int>(error_id) << \" -> \" << cudaGetErrorString(error_id) << std::endl;
          // std::cout << \"? \";
   //       exit(EXIT_FAILURE);
        }   
      }
    }
  ")
  execute_process(COMMAND "${CMAKE_CUDA_COMPILER}" "test.cu" WORKING_DIRECTORY ${tmpdir} RESULT_VARIABLE status ERROR_VARIABLE msg)
  if (NOT status EQUAL 0) 
    message(FATAL_ERROR "${pfx}: nvcc failed\n ${msg}")
  endif()
  execute_process(COMMAND "./a.out" WORKING_DIRECTORY ${tmpdir} RESULT_VARIABLE status OUTPUT_VARIABLE msg ERROR_VARIABLE err)
  message(STATUS "error output of CUDA devices detection: ${err}")
  if (NOT status EQUAL 0)
    message(FATAL_ERROR "${pfx}: test program failed")
  endif()
  message(STATUS "output of CUDA devices detection:\n ${msg}")
  # read device capabilites from the file created by test.cu
  file(READ "${tmpdir}/device_capability.txt" msg)
  separate_arguments(msg)
  list(REMOVE_ITEM msg "")
  list(LENGTH msg gpu_count)
  message(STATUS "number of CUDA devices detected: ${gpu_count}")
  if(gpu_count GREATER 1 OR LIBCLOUDPHXX_FORCE_MULTI_CUDA)
    # compile the multi_CUDA backend
    message(STATUS "Detected more than 1 GPU or LIBCLOUDPHXX_FORCE_MULTI_CUDA set, the multi_CUDA backend will be built.")
    target_sources(cloudphxx_lgrngn PRIVATE src/lib_multicuda.cu)
    target_compile_definitions(cloudphxx_lgrngn PRIVATE MULTICUDA_FOUND)
  else()
    message(STATUS "Number of GPUs detected <=1, the multi_CUDA backend will not be built. You can force compilation of this backend by setting LIBCLOUDPHXX_FORCE_MULTI_CUDA.")
  endif()
  list(REMOVE_DUPLICATES msg)

  if(msg)
    # find the lowest compute capability and set it as virtual architecture, see:
    # https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-compilation
    list(SORT msg)
    list(GET msg 0 VIRT)
#    if (VIRT STREQUAL "?")
#      set (VIRT "35")
#    endif()
    set_target_properties(cloudphxx_lgrngn PROPERTIES CUDA_ARCHITECTURES ${VIRT}-virtual)
#    set_property(TARGET  PROPERTY CUDA_ARCHITECTURES 70-real 72-virtual)

#    target_compile_options(cloudphxx_lgrngn PRIVATE  $<$<COMPILE_LANGUAGE:CUDA>: --gpu-architecture=compute_${VIRT}>)

    #set(GPU_CODE_OPTIONS  --gpu-code=)
    foreach(code ${msg})
      message(STATUS "CUDA capability: ${code}")
      if (code STREQUAL "?")
        set (code "35")
        message(STATUS "CUDA capability check failed, assuming a default of ${code}")
      endif()
      set(arch ${code})
      if (arch EQUAL "21") # see https://devtalk.nvidia.com/default/topic/606195/-solved-nvcc-fatal-unsupported-gpu-architecture-compute_21-/
        set(arch "20")
      endif()
      set_property(TARGET cloudphxx_lgrngn APPEND PROPERTY CUDA_ARCHITECTURES ${arch}-real)
      #set_target_properties(cloudphxx_lgrngn PROPERTIES CUDA_ARCHITECTURES ${arch}-real)
      #set(GPU_CODE_OPTIONS ${GPU_CODE_OPTIONS}sm_${arch},)
    endforeach(code)
#    target_compile_options(cloudphxx_lgrngn PRIVATE  $<$<COMPILE_LANGUAGE:CUDA>:${GPU_CODE_OPTIONS}>)
  else()
    message(STATUS "CUDA capability check didn't find any devices")
  endif()

  unset(pfx)
  unset(tmpdir)
  unset(msg)
  unset(status)

  # Release with debug info mode cuda flags, shouldn't cmake push CXX compiler flags as -Xcompiler=... ?
  if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    target_compile_options(cloudphxx_lgrngn PRIVATE  $<$<COMPILE_LANGUAGE:CUDA>: -O3 -use_fast_math -Xcompiler=-Ofast,-march=native,-fopenmp>)
  # Debug mode cuda flags
  elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(cloudphxx_lgrngn PRIVATE  $<$<COMPILE_LANGUAGE:CUDA>: -g -DTHRUST_DEBUG -lineinfo -Xcompiler=-fopenmp,-g,-Og,-rdynamic,-Wfatal-errors>)
  # Release mode cuda flags
  else()
    target_compile_options(cloudphxx_lgrngn PRIVATE  $<$<COMPILE_LANGUAGE:CUDA>: -DNDEBUG -O3 -use_fast_math -Xcompiler=-Ofast,-march=native,-DNDEBUG,-fopenmp>)
  endif()

  if(USE_MPI)
    target_compile_options(cloudphxx_lgrngn PRIVATE  $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler -DUSE_MPI>)
  endif()

  target_compile_definitions(cloudphxx_lgrngn PRIVATE CUDA_FOUND)

  # why was this added in the MPI PR? commented-out for now
  # set(CMAKE_REQUIRED_FLAGS " -fopenmp")

endif()

############################################################################################
# Thrust, location of Thrust can be hinted by setting THRUST_INCLUDE_DIR
find_package(Thrust REQUIRED)

#if (NOT THRUST_INCLUDE_DIR)
#  message(FATAL_ERROR "Thrust library not found.
#
#* To install Thrust, please try:
#*   Debian/Ubuntu: sudo apt-get install libthrust-dev
#*   Fedora: TODO
#*   Homebrew: TODO
#  ")
#endif()

# include thrust found here instead of the one shipped with cuda
#target_include_directories(cloudphxx_lgrngn PRIVATE ${THRUST_INCLUDE_DIR})

############################################################################################
# Boost libraries
find_package(Boost COMPONENTS) # we only need header-only boost libs
if(Boost_FOUND)
  target_link_libraries(cloudphxx_lgrngn
    PUBLIC 
      Boost::boost
  )
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
#message(STATUS "Testing if Boost ODEINT version >= 1.58")
#set(pfx "boost odeint check")
#execute_process(COMMAND "mktemp" "-d" "/tmp/tmp.XXX" RESULT_VARIABLE status OUTPUT_VARIABLE tmpdir)
#if (NOT status EQUAL 0)
#  message(FATAL_ERROR "${pfx}: mkdtemp failed")
#endif()
#file(WRITE "${tmpdir}/test.cpp" "
#  #define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CPP
#
#  #include <thrust/system/cpp/vector.h>
#
#  #include <boost/numeric/odeint.hpp>
#  #include <boost/numeric/odeint/external/thrust/thrust.hpp>
#
#  struct rhs
#  {
#    void operator()(
#      const thrust::cpp::vector<float> &psi,
#      thrust::cpp::vector<float> &dot_psi,
#      const float /* t */
#    )
#    {
#      assert(psi.size() == dot_psi.size());
#    }
#  };
#
#  int main()
#  {
#    boost::numeric::odeint::euler<
#      thrust::cpp::vector<float>, // state_type
#      float,                      // value_type
#      thrust::cpp::vector<float>, // deriv_type
#      float,                      // time_type
#      boost::numeric::odeint::thrust_algebra,
#      boost::numeric::odeint::thrust_operations
#    > chem_stepper;
#
#    thrust::cpp::vector<float> v(2);
#    chem_stepper.do_step(rhs(), v, 0, 1);
#  }
#")
#execute_process(COMMAND "${CMAKE_CXX_COMPILER}" "-std=c++14" "test.cpp" "-I${Boost_INCLUDE_DIR}" "-I${THRUST_INCLUDE_DIR}" WORKING_DIRECTORY ${tmpdir} RESULT_VARIABLE status ERROR_VARIABLE msg)
#if (NOT status EQUAL 0)
#  message(FATAL_ERROR "${pfx}: c++ compiler failed\n ${msg}")
#endif()
#execute_process(COMMAND "./a.out" WORKING_DIRECTORY ${tmpdir} RESULT_VARIABLE status OUTPUT_VARIABLE msg)
#if (NOT status EQUAL 0)
#  message(FATAL_ERROR "${pfx}: test program failed, install Boost odeint version >= 1.58")
#endif()
#unset(pfx)
#unset(tmpdir)
#unset(msg)
#unset(status)

# generate a header file with git revision id
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
  add_custom_target(git_revision.h ALL
   git log -1 "--format=format:#define LIBCLOUDPHXX_GIT_REVISION \"%H\"%n" HEAD > include/libcloudph++/git_revision.h
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} VERBATIM
  )
endif()

option(LIBCLOUDPHXX_DISABLE_BINDINGS "Disable compilation of Python bindings" OFF)

add_subdirectory(include) 

if(NOT LIBCLOUDPHXX_DISABLE_BINDINGS) # some tests require bindings, so we disable them too
  enable_testing()
  add_subdirectory(bindings)
  add_subdirectory(tests)
endif()

############################################################################################

# create a header-only library target (enough for blk_1m blk_2m and common)
add_library(cloudphxx_headers INTERFACE)

target_include_directories(cloudphxx_headers 
  INTERFACE
    $<INSTALL_INTERFACE:include>
)

target_compile_features(cloudphxx_headers INTERFACE cxx_std_14)

target_link_libraries(cloudphxx_headers
  INTERFACE
    Boost::boost
)

############################################################################################
# installation
set_target_properties(cloudphxx_lgrngn PROPERTIES DEBUG_POSTFIX _dbg RELWITHDEBINFO_POSTFIX _relwithdbg RELWITHDEBINFOPORTABLE_POSTFIX _relwithdbg)

# install files and export targets
install(TARGETS cloudphxx_lgrngn cloudphxx_headers EXPORT cloudphxx-targets 
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)

install(EXPORT cloudphxx-targets NAMESPACE clphxx:: DESTINATION share/libcloudph++)

#install(TARGETS cloudphxx_headers EXPORT cloudphxx_headers-targets)
#  INCLUDES DESTINATION include) # these files were already installed by lgrngn

#install(EXPORT cloudphxx_headers-targets NAMESPACE clphxx:: DESTINATION share/libcloudph++)

# generate and install a config file
include(CMakePackageConfigHelpers)

configure_package_config_file(
  libcloudph++-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/libcloudph++-config.cmake"
  INSTALL_DESTINATION 
    share/libcloudph++
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/libcloudph++-config.cmake"
  DESTINATION
    share/libcloudph++
)
