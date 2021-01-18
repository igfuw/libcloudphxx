#include "lib.hpp"

#include <thrust/system/cuda/vector.h>
namespace thrust_device = ::thrust::cuda;

#include "detail/gpu_assert.hpp"
#include "particles.tpp"

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
  };
};
