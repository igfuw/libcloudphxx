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
    void particles_t<real_t, device>::impl::post_adding_SD()
    {   
      thrust_device::vector<real_t> &drv = drv_gp->get();

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

      // update vt of new SD
      hskpng_vterm_invalid();

      // init _old values in per-particle substepping
      init_sstp();

      drv_gp.reset(); // release the tmp array that stored drv
    }
  };  
};
