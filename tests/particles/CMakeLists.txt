set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CMAKE_SOURCE_DIR}") # TODO: handle it correctly!!!
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -I${CMAKE_SOURCE_DIR}") # TODO: handle it correctly!!!

find_package(CUDA QUIET)
set(CUDA_PROPAGATE_HOST_FLAGS OFF) # -std=c++11 will not work

# trying to allow runtime choice betwee CUDA, CPP and OpenMP backends
set(files "")
set(files "${files};test_particles.cpp")
set(files "${files};particles_cpp.cpp")

if (OPENMP_FOUND)
  set(files "${files};particles_openmp.cpp")
endif()

if (CUDA_FOUND)
  set(files "${files};particles_cuda.cu")
  add_definitions(-DCUDA_FOUND)
  cuda_add_executable(test_particles ${files})
else()
  add_executable(test_particles ${files})
endif()
