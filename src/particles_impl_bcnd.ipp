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
      template <typename n_t, typename real_t>
      struct flag
      {
        BOOST_GPU_ENABLED
        n_t operator()(const real_t &)
        {
          return 0;
        }
      };

      struct count_vol
      {   
        template <typename n_t, typename real_t>
        BOOST_GPU_ENABLED
        real_t operator()(const n_t &n, const real_t &rw2)
        {
          return n * pow(rw2, real_t(3./2.));
        }
      };  
  
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
        case 1:
        {
          // hardcoded periodic boundary in x! (TODO - as an option)
          // when working on a single GPU simply apply bcond
          if(opts_init.dev_count < 2)
          {
            thrust::transform(
              x.begin(), x.end(),
              x.begin(),
              detail::periodic<real_t>(opts_init.x0, opts_init.x1)
            );
          }
          // more than one GPU - save ids of particles that need to be copied left/right
          else
          {
	    namespace arg = thrust::placeholders;
            // use i and k as temp storage - after bcond they are invalid anyway
            // multi_CUDA works only for 2D and 3D
            thrust_device::vector<thrust_size_t> &lft_id(i);
            thrust_device::vector<thrust_size_t> &rgt_id(k);

            // save ids of SDs to copy
            lft_count = thrust::copy_if(
              zero, zero+n_part,
              x.begin(),
              lft_id.begin(),
              arg::_1 < opts_init.x0
            ) - lft_id.begin();

            rgt_count = thrust::copy_if(
              zero, zero+n_part,
              x.begin(),
              rgt_id.begin(),
              arg::_1 >= opts_init.x1
            ) - rgt_id.begin();

            if(lft_count > in_n_bfr.size() || rgt_count > in_n_bfr.size())
              throw std::runtime_error("Overflow of the in/out buffer\n"); // TODO: resize buffers?
          }

          // hardcoded periodic boundary in y! (TODO - as an option)
          if (n_dims == 3)
          {
	    thrust::transform(
	      y.begin(), y.end(),
	      y.begin(),
	      detail::periodic<real_t>(opts_init.y0, opts_init.y1)
	    );
          }

          if (n_dims > 1)
          {
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

              thrust_device::vector<real_t> &n_filtered(tmp_device_real_part);

              thrust::fill(n_filtered.begin(), n_filtered.end(), 0.);

              thrust::transform_if(
                n.begin(), n.end(),               // input 1
                rw2.begin(),                      // input 2
                z.begin(),                        // stencil
                n_filtered.begin(),               // output
                detail::count_vol(),              // operation
                arg::_1 < opts_init.z0            // condition
              );
  
              ret = 4./3. * thrust::reduce(n_filtered.begin(), n_filtered.end())
#if !defined(__NVCC__)
              * pi<real_t>();
#else
              * CUDART_PI;
#endif

              // zero-out multiplicities
	      thrust::transform_if(   
		z.begin(), z.end(),          // input 
		n.begin(),                   // output
		detail::flag<n_t, real_t>(), // operation (zero-out)
		arg::_1 < opts_init.z0       // condition
	      );
	    }
          }
          break; 
        }
        case 0: break;
        default: assert(false);
      }

      return ret;
    }
  };  
};
