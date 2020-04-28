#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      enum bcond_t
      {
        sharedmem,    // copy to the same device
        distmem_cuda, // copy to another device on the same node
        distmem_mpi,  // copy to another device on another node
        open          // remove the SD
      };
    }; 
  };
};
