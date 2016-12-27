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
    void particles_t<real_t, device>::impl::init_dry_dry_sizes(
      real_t radius
    )
    {
      real_t rad3 = radius * radius * radius;
      thrust::fill(rd3.begin(), rd3.end(), rad3);
    }
  };
};
