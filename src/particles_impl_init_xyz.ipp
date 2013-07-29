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
    // init_xyz
    template <typename real_t, int device>
    void particles<real_t, device>::impl::init_xyz()
    {
      // TODO: wouldn't it be simpler to call a helper method 3 times?
      thrust_device::vector<real_t> 
                  *v[3] = {&x,      &y,      &z      };
      const int    n[3] = { opts.nx, opts.ny, opts.nz};
      const real_t d[3] = { opts.dx, opts.dy, opts.dz};
     
      for (int ix = 0; ix < 3; ++ix)
      {
        if (n[ix] == 0) continue;

        // memory allocation
        v[ix]->resize(n_part);

        // tossing random numbers [0,1] 
        rand_u01(n_part);

	// shifting from [0,1] to [0,nx*dx] 
        {
          using namespace thrust::placeholders;
	  thrust::transform(
	    u01.begin(), 
	    u01.end(), 
	    v[ix]->begin(), 
	    _1 * n[ix] * d[ix]
	  );
        }
      }
    }
  };
};
