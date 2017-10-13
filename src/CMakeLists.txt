find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# allowing runtime choice between CUDA, CPP and OpenMP backends
set(files "")
set(files "${files};lib.cpp")

set(files "${files};lib_cpp.cpp")
if (OPENMP_FOUND)
  set(files "${files};lib_omp.cpp")
endif()

if (CUDA_FOUND)
  set(files "${files};lib_cuda.cu")
  cuda_add_library(cloudphxx_lgrngn SHARED ${files})
  target_link_libraries(cloudphxx_lgrngn ${CUDA_curand_LIBRARY})
else()
  add_library(cloudphxx_lgrngn SHARED ${files})
endif()

set_target_properties(cloudphxx_lgrngn PROPERTIES DEBUG_POSTFIX _dbg RELWITHDEBINFO_POSTFIX _relwithdbg)
install(TARGETS cloudphxx_lgrngn LIBRARY DESTINATION lib)
