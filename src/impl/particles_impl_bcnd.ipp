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

      template<class real_t>
      struct count_vol
      {   
        real_t exponent;
        count_vol(real_t exponent) : exponent(exponent){};
        template <typename tuple>
        BOOST_GPU_ENABLED
        real_t operator()(const tuple &tup)
        {
          return 4./3. 
#if !defined(__NVCC__)
            * pi<real_t>()
#else
            * CUDART_PI
#endif
            * thrust::get<0>(tup)    // n
            * pow(
                thrust::get<1>(tup),  // radius at some power
                exponent);
        }
      };  

      template<class real_t>
      struct count_mass
      {   
        template <typename tuple>
        BOOST_GPU_ENABLED
        real_t operator()(const tuple &tup)
        {
          return thrust::get<0>(tup)  *  // n
                 thrust::get<1>(tup);    // chem_mass
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
    void particles_t<real_t, device>::impl::bcnd()
    {   
      switch (n_dims)
      {
        case 3:
        case 2:
        case 1:
        {
          // hardcoded periodic boundary in x! (TODO - as an option)
          // when working on a shared memory system, simply apply bcond
          if(!distmem())
          {
            thrust::transform(
              x.begin(), x.end(),
              x.begin(),
              detail::periodic<real_t>(opts_init.x0, opts_init.x1)
            );
          }
          // distributed memory - save ids of particles that need to be copied left/right
          else
          {
	          namespace arg = thrust::placeholders;

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

            const auto no_of_n_vctrs_copied(int(1));
            const auto no_of_real_vctrs_copied(distmem_real_vctrs.size());

            if(lft_count*no_of_n_vctrs_copied > in_n_bfr.size() || rgt_count*no_of_n_vctrs_copied  > in_n_bfr.size())
            {
              n_t new_size = lft_count > rgt_count ?
                               1.1 * lft_count : 
                               1.1 * rgt_count;

              std::cerr << "Overflow of the buffer, bfr size: " << in_n_bfr.size() << " to be copied left: " << lft_count << " right: " << rgt_count << "; resizing to: " << new_size << std::endl;

              in_n_bfr.resize(no_of_n_vctrs_copied * new_size);    
              out_n_bfr.resize(no_of_n_vctrs_copied * new_size);

              in_real_bfr.resize(no_of_real_vctrs_copied * new_size);
              out_real_bfr.resize(no_of_real_vctrs_copied * new_size);
            }
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

              // copy n of SDs that are out of the domain (otherwise remains n_filtered=0)
              thrust::transform_if(
                n.begin(), n.end(),               // input 1
                z.begin(),                        // stencil
                n_filtered.begin(),               // output
                thrust::identity<n_t>(),          // operation
                arg::_1 < opts_init.z0            // condition
              );

              // add total liquid water volume that fell out in this step
              output_puddle[outliq_vol] += 
                thrust::transform_reduce(
                  thrust::make_zip_iterator(thrust::make_tuple(
                    n_filtered.begin(), rw2.begin())),           // input start
                  thrust::make_zip_iterator(thrust::make_tuple(
                    n_filtered.begin(), rw2.begin())) + n_part,  // input end
                  detail::count_vol<real_t>(3./2.),              // operation
                  real_t(0),                                     // init val
                  thrust::plus<real_t>()
                );

              // add total dry volume that fell out in this step
              output_puddle[outdry_vol] += 
                thrust::transform_reduce(
                  thrust::make_zip_iterator(thrust::make_tuple(
                    n_filtered.begin(), rd3.begin())),           // input start
                  thrust::make_zip_iterator(thrust::make_tuple(
                    n_filtered.begin(), rd3.begin())) + n_part,  // input end
                  detail::count_vol<real_t>(1.),                 // operation
                  real_t(0),                                     // init val
                  thrust::plus<real_t>()
                );

              if(opts_init.chem_switch)
              {
                for (int i = 0; i < chem_all; ++i)
                  output_puddle[static_cast<output_t>(i)] += 
                    thrust::transform_reduce(
                      thrust::make_zip_iterator(thrust::make_tuple(
                        n_filtered.begin(), chem_bgn[i])),           // input start
                      thrust::make_zip_iterator(thrust::make_tuple(
                        n_filtered.end(), chem_end[i])),             // input end
                      detail::count_mass<real_t>(),                  // operation
                      real_t(0),                                     // init val
                      thrust::plus<real_t>()
                    );
              }

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
    }
  };  
};
