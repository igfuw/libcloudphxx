#if defined(CUDA_FOUND)
#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CUDA
#else
#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CPP
#endif

#include <thrust/system/omp/vector.h>
#define thrust_device ::thrust::omp // TODO: change to namespace thrust_device = ::thrust::omp?

#include "particles.tpp"

namespace libcloudphxx
{ 
  namespace common
  {
    namespace prtcls
    {
      // instantiation 
      template class particles<float, omp>;
      template class particles<double, omp>;
    };
  };
};
