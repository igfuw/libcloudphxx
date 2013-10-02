// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/theta_dry.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      /// @brief returns ret_t(x/c) 
      template <typename arg_t, typename ret_t>
      struct divide_by_constant_and_cast
      {
        arg_t c;
        divide_by_constant_and_cast(arg_t c) : c(c) {}

        __device__
        ret_t operator()(arg_t x) 
        { 
          return ret_t(x/c); 
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_ijk()
    {   
      // helper functor
      struct {
        void operator()(
          thrust_device::vector<real_t> &vx,
          thrust_device::vector<thrust_size_t> &vi,
          const real_t &vd
        ) {
	  thrust::transform(
	    vx.begin(), vx.end(),                                // input
	    vi.begin(),                                          // output
	    detail::divide_by_constant_and_cast<real_t, int>(vd) // operation
	  );
        }
      } helper;
      
      if (opts_init.nx != 0) helper(x, i, opts_init.dx);
      if (opts_init.ny != 0) helper(y, j, opts_init.dy);
      if (opts_init.nz != 0) helper(z, k, opts_init.dz);

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
            _1 * opts_init.nz + _2   // assuming z varies first (as in many other cases)
          );
          break;
        case 3:
          assert(false); // TODO!
          break;
        default:
          assert(false);
      }
      
      // flagging that particles are no longer sorted 
      sorted = false;
    }   
  };  
};
