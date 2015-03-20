#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      //return index corresponding to given radius, only works for files supplied by Jon
      template<class real_t, class n_t>
      BOOST_GPU_ENABLED
      n_t kernel_index(const real_t &R, const real_t &R_max = 1e10)
      {
        n_t ret;
#if !defined(__NVCC__)
        R <=100. ? ret = round(R) : ret = 100 + round((R-100.)/10.);
#else
        R <=100. ? ret = roundf(R) : ret = 100 + roundf((R-100.)/10.);
#endif
        if (R>R_max) ret = kernel_index<real_t, n_t> (R_max, R_max);
        return ret;
      }

      //return index of kernel_params corresponding to two given indices from Jon's files
      template<class n_t>
      BOOST_GPU_ENABLED
      n_t kernel_vector_index(const n_t &i, const n_t &j, const n_t &n_user_params = 0)
      {
        if(i >= j)
          return 0.5 * (i-1) * i + j + n_user_params;
        else
          return 0.5 * (j-1) * j + i + n_user_params;
      }
    }
  }
}
