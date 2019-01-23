#include <thrust/system/cpp/vector.h>
namespace thrust_device = ::thrust::cpp;

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
    template class particles_t<float, serial>;
    template class particles_t<double, serial>;
  };
};
