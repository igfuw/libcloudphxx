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
    void particles_t<real_t, device>::impl::rw_mom3_ante_change()
    {   
      // --- calc liquid water content before src ---
      hskpng_sort();
      reset_guardp(drw_mom3_gp, tmp_device_real_cell); 
      thrust_device::vector<real_t> &drw_mom3 = drw_mom3_gp->get();

      moms_all();
      moms_calc(rw2.begin(), real_t(3./2.));

      // drw_mom3 = -rw_mom3 ante change
      if(count_n!=n_cell)  
        thrust::fill(drw_mom3.begin(), drw_mom3.end(), real_t(0.));

      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drw_mom3.begin(), count_ijk.begin()), // output
        thrust::negate<real_t>()
      );
    }
  };  
};
