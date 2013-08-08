find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

find_package(CUDA QUIET)
set(CUDA_PROPAGATE_HOST_FLAGS OFF) # -std=c++11 will not work
# TODO: propagate optimisation flags to nvcc

# allowing runtime choice between CUDA, CPP and OpenMP backends
set(files "")
set(files "${files};lib.cpp")

set(files "${files};lib_cpp.cpp")
if (OPENMP_FOUND)
  set(files "${files};lib_omp.cpp")
endif()

# TODO: detecd CUDA arch (native functionality would be best?)
if (CUDA_FOUND)
  set(CUDA_NVCC_FLAGS "-I${CMAKE_SOURCE_DIR}/include ${CUDA_NVCC_FLAGS}")
  set(files "${files};lib_cuda.cu")
  add_definitions(-DCUDA_FOUND)
  cuda_add_library(cloudphxx_lgrngn SHARED ${files})
  target_link_libraries(cloudphxx_lgrngn ${CUDA_curand_LIBRARY})
else()
  add_library(cloudphxx_lgrngn SHARED ${files})
endif()

install(TARGETS cloudphxx_lgrngn LIBRARY DESTINATION lib)
