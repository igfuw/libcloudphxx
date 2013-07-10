#pragma once

#if (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP)
#  undef THRUST_DEVICE_SYSTEM
#  include <thrust/system/omp/vector.h>
namespace libcloudphxx { namespace common { namespace prtcls {
      namespace thrust_device = ::thrust::omp;
}}}
#elif (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA)
#  undef THRUST_DEVICE_SYSTEM
#  include <thrust/system/cuda/vector.h>
namespace libcloudphxx { namespace common { namespace prtcls {
      namespace thrust_device = ::thrust::cuda;
}}}
#elif (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CPP)
#  undef THRUST_DEVICE_SYSTEM
#  include <thrust/system/cpp/vector.h>
namespace libcloudphxx { namespace common { namespace prtcls {
      namespace thrust_device = ::thrust::cpp;
}}}
#else
#  error unknown or unspecified THRUST_DEVICE_SYSTEM
#endif

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      typedef thrust_device::vector<int>::size_type thrust_size_t;

//#if !defined(NDEBUG)
      namespace debug
      {
        template <typename real_t>
	void print(const thrust_device::vector<real_t> &v)
        {
	  thrust::copy(v.begin(), v.end(), std::ostream_iterator<real_t>(std::cerr, "\t"));
          std::cerr << std::endl;
        }
      };
//#endif
    };
  };
};
