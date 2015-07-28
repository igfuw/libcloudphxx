// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct periodic
      { 
        real_t a, b;

        periodic(real_t a, real_t b) : a(a), b(b) {}

        BOOST_GPU_ENABLED
        real_t operator()(real_t x)
        {
          return a + fmod((x-a) + (b-a), b-a); // this should call CUDA's fmod!
        }
      };
    };

    template <typename real_t, backend_t device>
    real_t particles_t<real_t, device>::impl::bcnd()
    {   
      real_t ret = 0;

      switch (n_dims)
      {
        case 3:
        case 2:
        {
          // hardcoded periodic boundary in x! (TODO - as an option)
          thrust::transform_if(
            x.begin(), x.end(),
            sd_stat.begin(),
            x.begin(),
            detail::periodic<real_t>(opts_init.x0, opts_init.x1),
            detail::is_active()
          );

          // hardcoded periodic boundary in y! (TODO - as an option)
          if (n_dims == 3)
          {
	    thrust::transform_if(
	      y.begin(), y.end(),
              sd_stat.begin(),
	      y.begin(),
	      detail::periodic<real_t>(opts_init.y0, opts_init.y1),
              detail::is_active()
	    );
          }

          // hardcoded "open" boudary at the top of the domain 
          // (just for numerical-error-sourced out-of-domain particles)
          {
            namespace arg = thrust::placeholders;
	    thrust::transform_if(
	      z.begin(), z.end(),          // input - arg
	      sd_stat.begin(),             // output
              detail::deactivate<real_t>(),      // operation (mark it as inactive, TODO: mark as to_rcyc after rcyc is fixed to work with sd_stat)
	      arg::_1 >= opts_init.z1      // condition (note: >= seems important as z==z1 would cause out-of-range ijk)
	    );
          }

          // precipitation on the bottom edge of the domain
          //// first: count the volume of particles below the domain
          // TODO! (using tranform_reduce?)
          //// second: zero-out multiplicities so they will be recycled
          {
            namespace arg = thrust::placeholders;
	    thrust::transform_if(   
	      z.begin(), z.end(),          // input 
	      sd_stat.begin(),             // output
	      detail::deactivate<real_t>(),      // operation (make it inactive)
	      arg::_1 < opts_init.z0       // condition
	    );
          }
          break; 
        }
        case 1:
          assert(false && "TODO");
          break;
        case 0: break;
        default: assert(false);
      }

      return ret;
    }
  };  
};
