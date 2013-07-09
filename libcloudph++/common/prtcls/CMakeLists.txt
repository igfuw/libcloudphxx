find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

find_package(CUDA QUIET)
set(CUDA_PROPAGATE_HOST_FLAGS OFF) # -std=c++11 will not work

# trying to allow runtime choice betwee CUDA, CPP and OpenMP backends
set(files "")
set(files "${files};particles.cpp")
set(files "${files};particles_cpp.cpp")

if (OPENMP_FOUND)
  set(files "${files};particles_omp.cpp")
endif()

if (CUDA_FOUND)
  set(files "${files};particles_cuda.cu")
  add_definitions(-DCUDA_FOUND)
  cuda_add_library(cloudphxx_common_prtcls SHARED ${files})
  target_link_libraries(cloudphxx_common_prtcls ${CUDA_curand_LIBRARY})
else()
  add_library(cloudphxx_common_prtcls SHARED ${files})
endif()

install(FILES particles.hpp DESTINATION include/libcloudph++/common/prtcls)
install(TARGETS cloudphxx_common_prtcls LIBRARY DESTINATION lib)
