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
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_xyz()
    {
      thrust_device::vector<real_t> 
                  *v[3] = { &x,           &y,           &z           };
      const int    n[3] = { opts_init.nx, opts_init.ny, opts_init.nz };
      const real_t a[3] = { opts_init.x0, opts_init.y0, opts_init.z0 };
      const real_t b[3] = { opts_init.x1, opts_init.y1, opts_init.z1 };
     
      for (int ix = 0; ix < 3; ++ix)
      {
        if (n[ix] == 0) continue;

        // memory allocation
        v[ix]->resize(n_part);

        // tossing random numbers [0,1] 
        rand_u01(n_part);

	// shifting from [0,1] to [x0,x1] 
        {
          using namespace thrust::placeholders;
	  thrust::transform(
	    u01.begin(), 
	    u01.end(), 
	    v[ix]->begin(), 
	    a[ix] + _1 * (b[ix] - a[ix])
	  );
        }
      }
    }
  };
};
