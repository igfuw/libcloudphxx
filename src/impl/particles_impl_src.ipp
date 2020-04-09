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
    void particles_t<real_t, device>::impl::src(const real_t &dt)
    {   
      // --- calc liquid water content before src ---
      hskpng_sort(); 
      thrust_device::vector<real_t> &drv(tmp_device_real_cell);
      thrust::fill(drv.begin(), drv.end(), real_t(0.));

      moms_all();
      moms_calc(rw2.begin(), real_t(3./2.));

      // drv = - tot_vol_bfr
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
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

      if(!opts_init.src_dry_distros.empty())
        src_dry_distros(dt);

      if(!opts_init.src_dry_sizes.empty())
        src_dry_sizes(dt);
 
      // --- after source particles are no longer sorted ---
      sorted = false;

      // --- calc liquid water content after src ---
      hskpng_sort(); 
      moms_all();
      moms_calc(rw2.begin(), real_t(3./2.));

      // drv = tot_vol_after -tot_vol_bfr + dry_vol_bfr
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // 2nd arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );

      // drv = tot_vol_after - dry_vol_after - tot_vol_bfr + dry_vol_bfr
/*
      moms_calc(rd3.begin(), 1);
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // 2nd arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::minus<real_t>()
      );
*/

      // update th and rv
      update_th_rv(drv);

      // update count_ijk and count_num
      hskpng_count();

      // store sstp_old
      sstp_save();
    }
  };  
};
