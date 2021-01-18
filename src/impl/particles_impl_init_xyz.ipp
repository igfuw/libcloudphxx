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
    };

    // reused in source
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_xyz()
    {
      // get i, j, k from ijk 
      switch(n_dims)
      {
        case 3:
          namespace arg = thrust::placeholders;
          // y
          thrust::transform(
            ijk.begin() + n_part_old, ijk.end(), // input - first arg
            j.begin() + n_part_old,        // output
            (arg::_1 / opts_init.nz) % (opts_init.ny) // z varies first
          );
          // z
          thrust::transform(
            ijk.begin() + n_part_old, ijk.end(), // input - first arg
            k.begin() + n_part_old,        // output
            arg::_1 % (opts_init.nz)   // z varies first
          );
          // x
          thrust::transform(
            ijk.begin() + n_part_old, ijk.end(), // input - first arg
            i.begin() + n_part_old,        // output
            arg::_1 / (opts_init.nz * opts_init.ny)    // z and y vary first
          );
          break;
        case 2:
          // z
          thrust::transform(
            ijk.begin() + n_part_old, ijk.end(), // input - first arg
            k.begin() + n_part_old,        // output
            arg::_1 % (opts_init.nz)   // z varies first
          );
          // x
          thrust::transform(
            ijk.begin() + n_part_old, ijk.end(), // input - first arg
            i.begin() + n_part_old,        // output
            arg::_1 / (opts_init.nz)
          );
          break;
        case 1:
          thrust::copy(ijk.begin() + n_part_old, ijk.end(), i.begin() + n_part_old); // only x
        case 0:
          break;
        default:
          assert(false);
          break;
      }

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

        // tossing random numbers [0,1] 
        rand_u01(n_part_to_init);

	// shifting from [0,1] to random position within respective cell 
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
    }
  };
};
