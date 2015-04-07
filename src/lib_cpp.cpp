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
    template <typename real_t, backend_t backend>
    void particles_t<real_t, backend>::impl::sanity_checks()
    {   
      detail::opts_init_sanity_checks(opts_init);
    }  

    // instantiation 
    template class particles_t<float, serial>;
    template class particles_t<double, serial>;
  };
};
