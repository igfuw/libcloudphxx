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

      template <typename n_t, typename real_t>
      struct flag
      {
        BOOST_GPU_ENABLED
        n_t operator()(const real_t &)
        {
          return 0;
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
          assert(false && "TODO");
          break; 
        case 2:
        {
          // hardcoded periodic boundary in x! (TODO - as an option)
          thrust::transform(
            x.begin(), x.end(),
            x.begin(),
            detail::periodic<real_t>(opts_init.x0, opts_init.x1)
          );

          // hardcoded "open" boudary at the top of the domain 
          // (just for numerical-error-sourced out-of-domain particles)
          {
            namespace arg = thrust::placeholders;
	    thrust::transform_if(
	      z.begin(), z.end(),          // input - arg
	      n.begin(),                   // output
              detail::flag<n_t, real_t>(), // operation (zero-out, so recycling will take care of it)
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
	      n.begin(),                   // output
	      detail::flag<n_t, real_t>(), // operation (zero-out)
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
