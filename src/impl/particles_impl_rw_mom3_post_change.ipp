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
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::rw_mom3_post_change()
    {   
      thrust_device::vector<real_t> &drw_mom3 = drw_mom3_gp->get();

      moms_all();
      moms_calc(rw2.begin(), real_t(3./2.));

      // drw_mom3 = -rw_mom3_ante + rw_mom3_post
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drw_mom3.begin(), count_ijk.begin()), // 2nd arg
        thrust::make_permutation_iterator(drw_mom3.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );
    }
  };  
};
