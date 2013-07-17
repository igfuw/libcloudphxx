#include "lib.hpp"

// workarounding Thrust bug #383: (Thanks to Jared for suggestion!)
#include <thrust/detail/minmax.h> 
#include <thrust/system/cpp/execution_policy.h>

#include <thrust/system/cpp/vector.h>
namespace thrust_device = ::thrust::cpp;

#include "particles.tpp"

namespace libcloudphxx
{ 
  namespace lgrngn
  {
    // instantiation 
    template class particles<float, cpp>;
    template class particles<double, cpp>;
  };
};
