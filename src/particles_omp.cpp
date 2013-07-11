#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP

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
