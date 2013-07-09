#pragma once

#include <thrust/device_vector.h>

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      typedef thrust::device_vector<int>::size_type thrust_size_t;

//#if !defined(NDEBUG)
      namespace debug
      {
        template <typename real_t>
	void print(const thrust::device_vector<real_t> &v)
        {
	  thrust::copy(v.begin(), v.end(), std::ostream_iterator<real_t>(std::cerr, "\t"));
          std::cerr << std::endl;
        }
      };
//#endif
    };
  };
};
