// workarounding Thrust bug #383: (Thanks to Jared for suggestion!)
#include <thrust/system/omp/execution_policy.h>

#include <thrust/system/omp/vector.h>
namespace thrust_device = ::thrust::omp;

#include "particles.tpp"

namespace libcloudphxx
{ 
  namespace lgrngn
  {
    // instantiation 
    template class particles<float, omp>;
    template class particles<double, omp>;
  };
};
