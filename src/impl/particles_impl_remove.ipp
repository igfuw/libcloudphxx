// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/reduce.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::remove_rng(
      const real_t &min, const real_t &max, 
      const typename thrust_device::vector<real_t>::iterator &vec_bgn
    )
    {
      namespace arg = thrust::placeholders;

      // zero-out multiplicities to mark for removal
      thrust::transform_if(   
        vec_bgn, vec_bgn + n_part,       // input 
        n.begin(),                       // output
        detail::flag<n_t, real_t>(),     // operation (zero-out)
        arg::_1 > min && arg::_1 <= max  // condition
      );
    }
  };  
};
