#pragma once

// to workaround nvcc unused vars warning
// http://stackoverflow.com/questions/777261/avoiding-unused-variables-warnings-when-using-assert-in-a-release-build
#define _unused(x) do { (void)sizeof(x); } while(0)

namespace libcloudphxx
{
  namespace lgrngn
  {
    typedef thrust_device::vector<int>::size_type thrust_size_t;

//#if !defined(NDEBUG) // TODO (CMake defaults)
    namespace debug
    {
      template <typename vec_t>
      void print(const vec_t &v)
      {
	thrust::copy(v.begin(), v.end(), std::ostream_iterator<typename vec_t::value_type>(std::cerr, "\t"));
	std::cerr << std::endl;
      }
    };
//#endif
  };
};
