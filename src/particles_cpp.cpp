//#if defined(CUDA_FOUND)
//#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CUDA
//#else
#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
//#endif

#include <thrust/system/cpp/vector.h>
#define thrust_device ::thrust::cpp // TODO: namespace

#include "particles.tpp"

namespace libcloudphxx
{ 
  namespace common
  {
    namespace prtcls
    {
      // instantiation 
      template class particles<float, cpp>;
      template class particles<double, cpp>;
    };
  };
};
