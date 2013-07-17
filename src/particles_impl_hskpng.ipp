// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include "detail/functors_device.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, int device>
    void particles<real_t, device>::impl::hskpng()
    {   
      // 0-D case
      if (n_dims == 0) return;
  
      // TODO: do it in a loop over simensions...
 
      // computing i
      thrust::transform(
        x.begin(), x.end(),
        i.begin(),
        detail::divide_by_constant_and_cast<real_t, int>(opts.dx)
      );
    }   
  };  
};
