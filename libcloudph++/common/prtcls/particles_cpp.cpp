#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CPP

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
