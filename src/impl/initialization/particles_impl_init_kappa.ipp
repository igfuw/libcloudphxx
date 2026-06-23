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
      const real_t kappa,     // kappa of the soluble part
      const real_t soluble_fraction   // volume fraction of the soluble part
    )
    {
      // assuming that rd3 = rd3_sol + rd3_insol is already filled
      // const real_t rd3_insol = rd_insol * rd_insol * rd_insol;

      thrust::fill(
        kpa.begin() + n_part_old, kpa.end(),
        kappa * soluble_fraction
      );
      
      // namespace arg = thrust::placeholders;
      // thrust::transform(
      //   rd3.begin() + n_part_old, rd3.end(),
      //   kpa.begin() + n_part_old,
      //   (arg::_1 - arg::_1 * soluble_fraction) / arg::_1 * kappa
      // );
    }
  };
};

