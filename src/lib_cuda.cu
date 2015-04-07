#include "lib.hpp"

// including it first not to require pthread option to nvcc
//#include <blitz/array.h>

// workarounding Thrust bug #383: (Thanks to Jared for suggestion!)
#include <thrust/system/cuda/execution_policy.h>

#include <thrust/system/cuda/vector.h>
namespace thrust_device = ::thrust::cuda;

#include "particles.tpp"

namespace libcloudphxx
{ 
  namespace lgrngn
  {
    template <typename real_t, backend_t backend>
    void particles_t<real_t, backend>::impl::sanity_checks()
    {   
      detail::opts_init_sanity_checks(opts_init);
    }  

    // instantiation 
    template class particles_t<float, CUDA>;
    template class particles_t<double, CUDA>; 
  };
};
