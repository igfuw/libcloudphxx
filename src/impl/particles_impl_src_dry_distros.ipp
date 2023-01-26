// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
// #include <limits>
#include <thrust/unique.h>
#include <thrust/binary_search.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // create new aerosol particles based on a size distribution
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::src_dry_distros(const real_t &dt)
    { 
      if(opts_init.src_type == src_t::matching)
        src_dry_distros_matching(dt);
      if(opts_init.src_type == src_t::simple)
        src_dry_distros_simple(dt);
    }
  };  
};
