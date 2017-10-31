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
    void particles_t<real_t, device>::impl::cond(
      const real_t &dt,
      const real_t &RH_max
    ) {   

      // --- calc liquid water content before cond ---
      hskpng_sort(); 
      thrust_device::vector<real_t> &drv(tmp_device_real_cell);

      // calculating the 3rd wet moment before condensation
      moms_all();
      moms_calc(rw2.begin(), real_t(3./2.));

      // permute-copying the result to -dm_3
      // fill with 0s if not all cells will be updated in the following transform
      if(count_n!=n_cell)  thrust::fill(drv.begin(), drv.end(), real_t(0.));
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::negate<real_t>()
      );

      // calculating drop growth in a timestep using backward Euler 
      thrust::transform(
        rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
        thrust::make_zip_iterator(      // input - 2nd arg
          thrust::make_tuple(
            thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
            thrust::make_permutation_iterator(rv.begin(), ijk.begin()),
            thrust::make_permutation_iterator(T.begin(), ijk.begin()),
            thrust::make_permutation_iterator(p.begin(), ijk.begin()),
            thrust::make_permutation_iterator(RH.begin(), ijk.begin()),
            thrust::make_permutation_iterator(eta.begin(), ijk.begin()),
            rd3.begin(),
            kpa.begin(),
            vt.begin()
          )
        ), 
        rw2.begin(),                    // output
        detail::advance_rw2<real_t>(dt, RH_max)
      );
      nancheck(rw2, "rw2 after condensation (no sub-steps");

      // calculating the 3rd wet moment after condensation
      moms_calc(rw2.begin(), real_t(3./2.));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) after condensation");

      // adding the third moment after condensation to dm_3
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // input - 2nd arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );

      // update th and rv according to changes in third specific wet moment
      update_th_rv(drv);
    }
  };  
};
