// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::remove_wet_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::remove_wet_rng, r_mi, r_mx);
    }
  };
};
