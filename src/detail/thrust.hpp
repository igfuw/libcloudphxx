#pragma once

// to workaround nvcc unused vars warning
// http://stackoverflow.com/questions/777261/avoiding-unused-variables-warnings-when-using-assert-in-a-release-build
#define _unused(x) do { (void)sizeof(x); } while(0)

namespace libcloudphxx
{
  namespace lgrngn
  {
    typedef thrust_device::vector<int>::size_type thrust_size_t;
  };
};
