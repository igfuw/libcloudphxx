#include "tpl_rw_t_wrapper.hpp"

#if defined(__NVCC__)
#  include <math_constants.h>
#endif

namespace libcloudphxx
{
  namespace lgrngn
  {
    using detail::tpl_rw_t_wrap;

    template <typename real_t, typename n_t>
    struct kernel_base
    {
      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_rw_t_wrap<real_t,n_t> &) const {return 0;}
    };

    //Golovin kernel
    template <typename real_t, typename n_t>
    struct kernel_golovin : kernel_base<real_t, n_t>
    {
      BOOST_GPU_ENABLED
      real_t calc(const tpl_rw_t_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };
#if !defined(__NVCC__)
        using std::abs;
        using std::pow;
        using std::max;
#endif
        real_t res =
#if !defined(__NVCC__)
        pi<real_t>()
#else
        CUDART_PI
#endif
        * real_t(2000.) //= 4/3 * (b=1500) to compare with Shima 2009, TODO: accept parameter from opts_init
        * max(
            thrust::get<n_a_ix>(tpl_wrap()),
            thrust::get<n_b_ix>(tpl_wrap())
          )
        * (
            pow(thrust::get<rw2_a_ix>(tpl_wrap()),real_t(3./2.)) +
            pow(thrust::get<rw2_b_ix>(tpl_wrap()),real_t(3./2.))
          );
        return res;
      }
    };

    //geometric kernel
    template <typename real_t, typename n_t>
    struct kernel_geometric : kernel_base<real_t, n_t>
    {
      BOOST_GPU_ENABLED
      real_t calc(const tpl_rw_t_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };
#if !defined(__NVCC__)
        using std::abs;
        using std::pow;
        using std::max;
#endif
        real_t res =
#if !defined(__NVCC__)
        pi<real_t>()
#else
        CUDART_PI
#endif
        * max(
            thrust::get<n_a_ix>(tpl_wrap()),
            thrust::get<n_b_ix>(tpl_wrap())
          )
        * abs(
            thrust::get<vt_a_ix>(tpl_wrap()) -
            thrust::get<vt_b_ix>(tpl_wrap())
          )
        * pow(
            sqrt(thrust::get<rw2_a_ix>(tpl_wrap())) +
            sqrt(thrust::get<rw2_b_ix>(tpl_wrap())),
            real_t(2)
          );
        return res;
      }
    };
  };
};
                            
