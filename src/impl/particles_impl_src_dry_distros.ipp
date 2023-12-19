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
    void particles_t<real_t, device>::impl::src_dry_distros(const real_t &dt, const dry_distros_t<real_t> &sdd)
    { 
      if (sdd.size() > 1)
        throw std::runtime_error("libcloudph++: src_dry_distros can only have a single kappa value.");

      if (opts_init.src_type == src_t::matching && !sdd.empty() &&
          sdd.begin()->first != opts_init.dry_distros.begin()->first) throw std::runtime_error("libcloudph++: For 'matching' CCN source, kappa of the source has to be the same as that of the initial profile (no kappa matching done)");

      if(opts_init.src_type == src_t::matching)
        src_dry_distros_matching(dt, sdd);
      if(opts_init.src_type == src_t::simple)
        src_dry_distros_simple(dt, sdd);
    }
  };  
};
