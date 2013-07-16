//#if defined(CUDA_FOUND)
//#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CUDA
//#else
#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
//#endif

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
