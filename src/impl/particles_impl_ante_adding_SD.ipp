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
    void particles_t<real_t, device>::impl::ante_adding_SD()
    {   
      // --- calc liquid water content before src ---
      hskpng_sort(); 
      thrust_device::vector<real_t> &drv(tmp_device_real_cell1); // NOTE: this can't be changed by any function called before a call to after_adding_SD...
      thrust::fill(drv.begin(), drv.end(), real_t(0.));

      moms_all();
      moms_calc(rw2.begin(), real_t(3./2.));

      // drv = - tot_vol_bfr
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n.get(),                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::negate<real_t>()
      );

      // drv = -tot_vol_bfr + dry_vol_bfr
/*
      moms_calc(rd3.begin(), 1);
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // 2nd arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );
*/
    }
  };  
};
