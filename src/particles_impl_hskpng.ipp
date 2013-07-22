// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include "detail/functors_device.hpp"
//#include "../include/common/theta_dry.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, int device>
    void particles<real_t, device>::impl::hskpng_ijk()
    {   
      // TODO: check if it works with nvcc and remove the old code
/*
      thrust_device::vector<real_t> 
                  *vx[3] = {&x,      &y,      &z      };
      thrust_device::vector<int> 
                  *vi[3] = {&i,      &j,      &k      };  
      const int    vn[3] = { opts.nx, opts.ny, opts.nz};
      const real_t vd[3] = { opts.dx, opts.dy, opts.dz};
*/
    
      // helper functor
      struct {
        void operator()(
          thrust_device::vector<real_t> &vx,
          thrust_device::vector<int> &vi,
          const real_t &vd
        ) {
	  thrust::transform(
	    vx.begin(), vx.end(),                          // input
	    vi.begin(),                                         // output
	    detail::divide_by_constant_and_cast<real_t, int>(vd) // operation
	  );
        }
      } helper;
      
      if (opts.nx != 0) helper(x, i, opts.dx);
      if (opts.ny != 0) helper(y, j, opts.dy);
      if (opts.nz != 0) helper(z, k, opts.dz);

/*
      // computing i, j & k into ijk
      for (int ix = 0; ix < 3; ++ix)
      {   
        if (vn[ix] == 0) continue;
 
        thrust::transform(
          vx[ix]->begin(), vx[ix]->end(),                          // input
          vi[ix]->begin(),                                         // output
          detail::divide_by_constant_and_cast<real_t, int>(vd[ix]) // operation
        );
      }
*/

      // raveling i, j & k into ijk
      switch (n_dims)
      {
        case 0: 
          ijk[0] = 0;
          break;
        case 1:
          thrust::copy(k.begin(), k.end(), ijk.begin());
          break;
        case 2:
          using namespace thrust::placeholders;
          thrust::transform(
            i.begin(), i.end(), // input - first arg
            k.begin(),          // input - second arg
            ijk.begin(),        // output
            _1 * opts.nz + _2   // assuming z varies first (as in many other cases)
          );
          break;
        case 3:
          assert(false); // TODO!
          break;
        default:
          assert(false);
      }
    }   

    template <typename real_t, int device>
    void particles<real_t, device>::impl::hskpng_Tpr()
    {   
      using namespace thrust::placeholders;
      // r  = rhod_rv / rhod;
      thrust::transform(rhod_rv.begin(), rhod_rv.end(), rhod.begin(), r.begin(), _1 / _2);
      // T  = common::theta_dry::T<real_t>(rhod_th, rhod);
      //thrust::transform(rhod_th.begin(), rhod_th.end(), rhod.begin(), T.begin(), common::theta_dry::T<real_t>(_1, _2));
      // p  = common::theta_dry::p<real_t>(rhod, r, T); 

    }
  };  
};
