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
    // save number of SDs in cells
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_count_num(const real_t &ratio)
    {
      if(n_dims > 0)
      {
        namespace arg = thrust::placeholders;
        // some cells may be used only partially in thr super-droplet method
        // e.g. when Lagrangian domain (x0, x1, etc...) is smaller than the 
        // Eulerian domain (0, nx*dx, etc...)
        // sd_conc defines number of SDs per Eulerian cell
        thrust::transform(dv.begin(), dv.end(), count_num.begin(), (ratio * real_t(opts_init.sd_conc) * arg::_1 / (opts_init.dx * opts_init.dy * opts_init.dz) + real_t(0.5))); 
      }
      // parcel setup
      else
        thrust::fill(count_num.begin(), count_num.end(), ratio * opts_init.sd_conc);
    }
  };
};
