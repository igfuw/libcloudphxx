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
    // init_xyz, to get uniform distribution in each cell
    // first n_cell SDs are distributed one per each cell,
    // then same with second n_cell particles, etc.
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_xyz()
    {
      thrust_device::vector<real_t> 
                  *v[3] = { &x,           &y,           &z           };
      const int    n[3] = { opts_init.nx, opts_init.ny, opts_init.nz };
      const real_t a[3] = { opts_init.x0, opts_init.y0, opts_init.z0 };
      const real_t d[3] = { opts_init.dx, opts_init.dy, opts_init.dz };
      thrust_device::vector<thrust_size_t> 
                  *ii[3] = { &i,           &j,           &k           };

      if(opts_init.sd_conc_mean > 0)
        thrust::fill(count_num.begin(), count_num.end(), opts_init.sd_conc_mean); // if using const_multi, count_num is already filled
 
      thrust::sequence(sorted_ijk.begin(), sorted_ijk.end()); // temporarily use sorted_ijk
      // get sorted_ijk = {0, 1,..., n_cell, 0, 1, ..., n_cell, ...}
      thrust::transform(sorted_ijk.begin(), sorted_ijk.end(), thrust::make_constant_iterator<thrust_size_t>(n_cell), sorted_ijk.begin(), thrust::modulus<thrust_size_t>());

      // get random keys
      rand_u01(n_part);
      {
        namespace arg = thrust::placeholders;
        // increment random keys by 1 after each n_cell keys
        thrust::transform(u01.begin(), u01.end(), thrust::make_counting_iterator(0), u01.begin(), arg::_1 + arg::_2 / n_cell);
      }

      // obtain a sequence of n_part/n_cell random permutations of the range [0, .., n_cell-1]
      thrust::sort_by_key(u01.begin(), u01.end(), sorted_ijk.begin());
     
      // get i, j, k from sorted_ijk 
      // i, j, k will be temporarily used
      switch(n_dims)
      {
        case(0):
          break;
        case(1):
          // z
          thrust::copy(sorted_ijk.begin(), sorted_ijk.end(), k.begin());
          break;
        case 2:
          namespace arg = thrust::placeholders;
          // x
          thrust::transform(
            sorted_ijk.begin(), sorted_ijk.end(), // input - first arg
            i.begin(),        // output
            arg::_1 / opts_init.nz   // assuming z varies first
          );
          // z
          thrust::transform(
            sorted_ijk.begin(), sorted_ijk.end(), // input - first arg
            k.begin(),        // output
            arg::_1 % opts_init.nz   // assuming z varies first
          );
          break;
        case 3:
          namespace arg = thrust::placeholders;
          // y
          thrust::transform(
            sorted_ijk.begin(), sorted_ijk.end(), // input - first arg
            j.begin(),        // output
            arg::_1 / (opts_init.nz * opts_init.nx)   // assuming z and x vary first
          );
          // x
          thrust::transform(
            sorted_ijk.begin(), sorted_ijk.end(), // input - first arg
            i.begin(),        // output
            arg::_1 % (opts_init.nz * opts_init.nx) / (opts_init.nz)   // assuming z varies first
          );
          // z
          thrust::transform(
            sorted_ijk.begin(), sorted_ijk.end(), // input - first arg
            k.begin(),        // output
            arg::_1 % (opts_init.nz * opts_init.nx) % (opts_init.nz)   // assuming z varies first
          );
          break;
        default:
          assert(false);
      }
      for (int ix = 0; ix < 3; ++ix)
      {
        if (n[ix] == 0) continue;

        // memory allocation
        v[ix]->resize(n_part);

        // tossing random numbers [0,1] 
        rand_u01(n_part);

	// shifting from [0,1] to random position within respective cell 
        {
          namespace arg = thrust::placeholders;
	  thrust::transform(
	    u01.begin(), 
	    u01.end(),
            ii[ix]->begin(), 
	    v[ix]->begin(), 
	    a[ix] + (arg::_1 + arg::_2) *d[ix]
	  );
        }
      }
    }
  };
};
