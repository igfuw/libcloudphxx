#pragma once

#include <thrust/device_vector.h>

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      typedef thrust::device_vector<int>::size_type thrust_size_t;
    };
  };
};
