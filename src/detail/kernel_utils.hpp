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
      n_t kernel_index(const real_t &R)
      {
        if( R<= 100.)
          return R;
        else
          return (100 + (R-100.) / 10.);
      }

      //return index of kernel_params corresponding to two given indices from Jon's files
      template<class n_t>
      BOOST_GPU_ENABLED
      n_t kernel_vector_index(const n_t &i, const n_t &j, const n_t &n_user_params = 0)
      {
        if(i>=j)
          return 0.5*i*(i+1) + j + n_user_params;
        else
          return 0.5*j*(j+1) + i + n_user_params;
      }
    }
  }
}
