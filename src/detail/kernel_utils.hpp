#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      //return index corresponding to given radius
      template<class n_t>
      BOOST_GPU_ENABLED
      int kernel_index(const n_t &R)
      {
        if( R<= 100.)
          return R;
        else
          return 100 + (R-100.) / 10.;
      }

      //return index of kernel_params corresponding to two given indices from Jon's files
      template<class n_t>
      BOOST_GPU_ENABLED
      thrust_size_t kernel_vector_index(const int &i, const int &j, const n_t &n_user_params = 0)
      {
        if(i>=j)
          return 0.5*i*(i+1) + j + n_user_params;
        else
          return 0.5*j*(j+1) + i + n_user_params;
      }
    }
  }
}
