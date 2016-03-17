#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      template <typename real_t>
      struct remote
      {   
        real_t lcl, rmt;

        remote(real_t lcl, real_t rmt) : lcl(lcl), rmt(rmt) {}

        BOOST_GPU_ENABLED
        real_t operator()(real_t x)
        {
          real_t res = rmt + x - lcl;
          if(res == rmt) res = nextafter(res, real_t(0.)); // in single precision, we used to get x=x1
          return res;
        }
      };
    }; 
  };
};
