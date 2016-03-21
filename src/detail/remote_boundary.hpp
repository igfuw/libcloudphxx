#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      enum bcond_t
      {
        sharedmem, distmem_cuda, distmem_mpi
      };

      template <typename real_t>
      struct remote
      {   
        real_t lcl, rmt;

        remote(real_t lcl, real_t rmt) : lcl(lcl), rmt(rmt) {}

        BOOST_GPU_ENABLED
        real_t operator()(real_t x)
        {
          real_t res = rmt + x - lcl;
          if(res == rmt) res = nextafter(res, real_t(0.)); // in single precision, we used to get x=x1, TODO: does it call CUDA's nextafter when used by multi_CUDA?
          return res;
        }
      };
    }; 
  };
};
