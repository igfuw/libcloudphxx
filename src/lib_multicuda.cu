#include "lib.hpp"

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

    template class particles_t<float, multi_CUDA>;
    template class particles_t<double, multi_CUDA>; 
  };
};
