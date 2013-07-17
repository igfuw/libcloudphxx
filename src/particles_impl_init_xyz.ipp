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
      thrust_device::vector<real_t> 
                  *v[3] = {&x,      &y,      &z      };
      const int    n[3] = { opts.nx, opts.ny, opts.nz};
      const real_t d[3] = { opts.dx, opts.dy, opts.dz};
     
      for (int i = 0; i < 3; ++i)
      {
        // memory allocation
        if (n[i] != 0) v[i]->resize(n_part);

        // tossing random numbers [0,1] 
        urand(n_part);

	// shifting from [0,1] to [0,nx*dx] 
        {
          using namespace thrust::placeholders;
	  thrust::transform(
	    u01.begin(), 
	    u01.end(), 
	    v[i]->begin(), 
	    _1 * n[i] * d[i]
	  );
        }
      }
    }
  };
};
