#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
//#if !defined(NDEBUG) // TODO (CMake defaults)
    namespace debug
    {
      template <typename vec_t>
      void print(const vec_t &v)
      {
        thrust::copy(v.begin(), v.end(), std::ostream_iterator<typename vec_t::value_type>(std::cerr, "\t"));
        std::cerr << std::endl;
      }

      template <typename vec_bgn_t, typename vec_end_t>
      void print(const vec_bgn_t &bgn, const vec_end_t &end)
      {
        thrust::copy(bgn, end, std::ostream_iterator<typename vec_bgn_t::value_type>(std::cerr, "\t"));
        std::cerr << std::endl;
      }

      template<class T>
      void print(ref_grid<T> v)
      {
        std::cerr << "normal grid values:" << std::endl;
        print(v.begin(), v.end());
        std::cerr << "refined grid values:" << std::endl;
        print(v.get_ref());
      }
    };
//#endif
  };
};
