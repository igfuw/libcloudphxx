//#if defined(CUDA_FOUND)
//#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CUDA
//#else
#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
//#endif

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
