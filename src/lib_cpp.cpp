#include "lib.hpp"

// workarounding Thrust bug #383: (Thanks to Jared for suggestion!)
#include <thrust/system/cpp/execution_policy.h>

#include <thrust/system/cpp/vector.h>
namespace thrust_device = ::thrust::cpp;

#include "particles.tpp"

namespace libcloudphxx
{ 
  namespace lgrngn
  {
    template <typename real_t, int backend>
    void particles_t<real_t, backend>::impl::sanity_checks()
    {   
    }  

    // instantiation 
    template class particles_t<float, cpp>;
    template class particles_t<double, cpp>;
  };
};
