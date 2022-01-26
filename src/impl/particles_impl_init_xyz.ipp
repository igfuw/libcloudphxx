// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <thrust/sequence.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<typename real_t>
      struct pos_lgrngn_domain
      // get a random position within ii-th cell taking into account Lagrangian domain
      {
        real_t p0, // lower bound of the Lagrangian domain
               p1, // upper bound of the Lagrangian domain
               dp; // cell size of the Eulerian grid

        pos_lgrngn_domain(real_t p0, real_t p1, real_t dp): p0(p0), p1(p1), dp(dp) {}

        BOOST_GPU_ENABLED
        real_t operator()(real_t u01, thrust_size_t ii) // random number [0,1), cell index in respective dimension
        {
#if !defined(__NVCC__)
          using std::min;
          using std::max;
#endif
        	
          return u01 * min(p1, (ii+1) * dp) + (1. - u01) * max(p0, ii * dp); 
        }
      };

     // returns cos(2*PI*arg)
     template <class real_t>
     struct cos2Pi : thrust::unary_function<const real_t&, real_t>
     {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &arg)
        {
#if !defined(__NVCC__)
          using std::cos;
#endif
        	
          return cos(2*arg*
#if !defined(__NVCC__)
            pi<real_t>()
#else
            CUDART_PI
#endif
          );
        }
      };

     // returns sin(2*PI*arg)
     template <class real_t>
     struct sin2Pi : thrust::unary_function<const real_t&, real_t>
     {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &arg)
        {
#if !defined(__NVCC__)
          using std::sin;
#endif
        	
          return sin(2*arg*
#if !defined(__NVCC__)
            pi<real_t>()
#else
            CUDART_PI
#endif
          );
        }
      };

      // returns arg_1 * sqrt(1 - arg_2^2)
      template <class real_t>
      struct mul_by_sqrt_1_min_square
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &arg1, const real_t &arg2)
        {
#if !defined(__NVCC__)
          using std::sqrt;
#endif
          return arg1 * sqrt(real_t(1) - arg2 * arg2);
        }
      };
    };

    // reused in source
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_xyz()
    {
      thrust_device::vector<real_t> 
                  *v[3] = { &x,           &y,           &z           };
      const int    n[3] = { opts_init.nx, opts_init.ny, opts_init.nz };
      const real_t a[3] = { opts_init.x0, opts_init.y0, opts_init.z0 };
      const real_t b[3] = { opts_init.x1, opts_init.y1, opts_init.z1 };
      const real_t d[3] = { opts_init.dx, opts_init.dy, opts_init.dz };
      thrust_device::vector<thrust_size_t> 
                  *ii[3] = { &i,           &j,           &k           };

      for (int ix = 0; ix < 3; ++ix)
      {
        if (n[ix] == 0) continue;

        // tossing random numbers
        rand_u01(n_part_to_init);

	// shifting from [0,1] to random position within respective cell 
  // TODO: now the rand range is [0,1), include this here
        {
          namespace arg = thrust::placeholders;
	  thrust::transform(
	    u01.begin(), 
	    u01.begin() + n_part_to_init,
            ii[ix]->begin() + n_part_old, 
	    v[ix]->begin() + n_part_old, 
            detail::pos_lgrngn_domain<real_t>(a[ix], b[ix], d[ix])
	  );
        }
      }

      // handle the initial pair separation option 
      // position of every second droplet seprated randomly by init_pair_separation distance [m] from the previous one
      // works only in 3D
      // NOTE: the second droplet of a pair can end up in a different cell!
      if(opts_init.init_pair_separation >= 0)
      {
        // random unit vector (e0,e1,e2)
        // h  [-1:1]
        // th [0:2*Pi]
        // e0 = sqrt(1. - h*h) * cos(th);
        // e1 = sqrt(1. - h*h) * sin(th);
        // e2 = h;
        
        thrust_device::vector<real_t> e0(tmp_device_real_part1);
        thrust_device::vector<real_t> e1(n_part_to_init); // TODO: use tmp_device_real_part2 (first ensure it is resized if init_pair_separation==True)
        thrust_device::vector<real_t> e2(tmp_device_real_part); // aka u01

        // fill u01 with random numbers [0,1] (angle th)
        rand_u01(n_part_to_init);

        namespace arg = thrust::placeholders;

        thrust::transform(
          u01.begin(),
          u01.begin() + n_part_to_init,
          e0.begin(),
          detail::cos2Pi<real_t>()
        );

        thrust::transform(
          u01.begin(),
          u01.begin() + n_part_to_init,
          e1.begin(),
          detail::sin2Pi<real_t>()
        );

        // fill e2 with random numbers [-1,1]
        rand_u01(n_part_to_init);
        thrust::transform(
          u01.begin(),
          u01.begin() + n_part_to_init,
          e2.begin(),
          real_t(2) * arg::_1 - real_t(1)
        );

        // multiply e0 and e1 by sqrt(1 - e2*e2)
        thrust_device::vector<real_t> *e[3] = {&e0, &e1, &e2};
        for(int i=0; i<2; ++i)
        {
          thrust::transform(
            e[i]->begin(),
            e[i]->begin() + n_part_to_init,
            e2.begin(),
            e[i]->begin(), // in place
            detail::mul_by_sqrt_1_min_square<real_t>()
          );
        }

        // set position of every second to pos_first + e * pair_separation
        for(int i=0; i<3; ++i)
        {
          thrust::transform_if(
            v[i]->begin(),                                       // arg 1 - position
            v[i]->begin() + n_part_to_init - 1,     
            e[i]->begin(),                                       // arg 2 - unit vector
            thrust::make_counting_iterator<thrust_size_t>(0),    // stencil
            v[i]->begin() + 1,                                   // output
            arg::_1 + arg::_2 * opts_init.init_pair_separation,  // operation
            arg::_1 % thrust_size_t(2) == thrust_size_t(0)       // predicate
          );
        }

        // handle potential out-of-bounds droplets
        bcnd();
      }
    }
  };
};
