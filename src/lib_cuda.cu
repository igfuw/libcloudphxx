#include "lib.hpp"

// including it first not to require pthread option to nvcc
//#include <blitz/array.h>

// workarounding Thrust bug #383: (Thanks to Jared for suggestion!)
#include <thrust/system/cuda/execution_policy.h>

#include <thrust/system/cuda/vector.h>
namespace thrust_device = ::thrust::cuda;

#include "particles.tpp"
#include "particles_multi_gpu.tpp"

namespace libcloudphxx
{ 
  namespace lgrngn
  {
    template <typename real_t, backend_t backend>
    void particles_t<real_t, backend>::impl::sanity_checks()
    {   
    }  

    // instantiation 
    template class particles_t<float, CUDA>;
    template class particles_t<double, CUDA>; 

    // TODO: move these to other file added if cmake detects more than 1 GPU?
    template class particles_t<float, multi_CUDA>;
    template class particles_t<double, multi_CUDA>; 
  };
};
