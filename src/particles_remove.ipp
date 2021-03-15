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
    // removes particles with (r_w >= r_min && r_w < r_max)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::remove_wet_rng(const real_t &r_min, const real_t &r_max)
    {
#if !defined(__NVCC__)
      using std::pow;
#endif
      pimpl->remove_rng(pow(r_min, 2), pow(r_max, 2), pimpl->rw2.begin());
    }
  };
};
