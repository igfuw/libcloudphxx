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
    // init number of SDs to be initialized per cell
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_count_num(const real_t &ratio)
    {
      thrust::fill(count_num.begin(), count_num.end(), ratio * opts_init.sd_conc);
    }
  };
};
