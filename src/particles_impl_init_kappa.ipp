// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_kappa(
      const real_t &kappa
    )
    {
      // filling kappas
      thrust::fill(kpa.begin() + n_part_old, kpa.end(), kappa);
    }
  };
};

