// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/theta_dry.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      /// @brief returns ret_t(x/c) 
      template <typename arg_t, typename ret_t>
      struct divide_by_constant_and_cast
      {
        arg_t c;
        divide_by_constant_and_cast(arg_t c) : c(c) {}

        BOOST_GPU_ENABLED
        ret_t operator()(arg_t x) 
        { 
          return ret_t(x/c); 
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_ijk()
    {   
      // helper functor
      struct {
        void operator()(
          thrust_device::vector<real_t> &vx,
          thrust_device::vector<thrust_size_t> &vi,
          thrust_device::vector<detail::sd_stat_t> &sd_stat,
          const real_t &vd
        ) {
	  thrust::transform_if(
	    vx.begin(), vx.end(),                                 // input
            sd_stat.begin(),                                      // stencil
	    vi.begin(),                                           // output
	    detail::divide_by_constant_and_cast<real_t, int>(vd), // operation
            detail::is_active_or_to_init()                           // predicate
	  );
        }
      } helper;
      
      if (opts_init.nx != 0) helper(x, i, sd_stat, opts_init.dx);
      if (opts_init.ny != 0) helper(y, j, sd_stat, opts_init.dy);
      if (opts_init.nz != 0) helper(z, k, sd_stat, opts_init.dz);

      // raveling i, j & k into ijk
      switch (n_dims)
      {
        case 0: 
          break;
        case 1:
          thrust::copy(k.begin(), k.end(), ijk.begin());
          break;
        case 2:
          namespace arg = thrust::placeholders;
          thrust::transform_if(
            i.begin(), i.end(), // input - first arg
            k.begin(),          // input - second arg
            sd_stat.begin(),
            ijk.begin(),        // output
            arg::_1 * opts_init.nz + arg::_2,   // assuming z varies first
            detail::is_active_or_to_init()
          );
          break;
        case 3:
          namespace arg = thrust::placeholders;
          thrust::transform_if(
            i.begin(), i.end(), // input - first arg
            j.begin(),          // input - second arg
            sd_stat.begin(),
            ijk.begin(),        // output
            arg::_1 * (opts_init.nz * opts_init.ny) + 
            arg::_2 * opts_init.nz,
            detail::is_active_or_to_init()
          );
          thrust::transform_if(
            ijk.begin(), ijk.end(),
            k.begin(),
            sd_stat.begin(),
            ijk.begin(), // in-place!
            arg::_1 + arg::_2,
            detail::is_active_or_to_init()
          );
          break;
        default:
          assert(false);
      }
      
      // flagging that particles are no longer sorted 
      sorted = false;
    }   
  };  
};
