#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      enum bcond_t
      {
        sharedmem,         // single process (no MPI) with serial / OpenMP / CUDA / multi_CUDA with 1 device 
        distmem_cuda_intr, // boundary between CUDA devices managed by one multi_CUDA process, except for the cyclic boundary on the edges
        distmem_cuda_extr, // a cyclic boundary on the edges between CUDA devices managed by one multi_CUDA process
        distmem_mpi        // boundary between different processes
      };

      inline bool bcond_is_distmem_cuda(bcond_t bcond) {return bcond == distmem_cuda_intr || bcond == distmem_cuda_extr;}
    }; 
  };
};
