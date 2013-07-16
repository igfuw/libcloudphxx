#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    typedef thrust_device::vector<int>::size_type thrust_size_t;

//#if !defined(NDEBUG) // TODO (CMake defaults)
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
