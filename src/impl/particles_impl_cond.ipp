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

    namespace detail
    {
      // Functor to select moment of only liquid droplets
      template<class real_t>
      class select_liq
      {
      public:
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, int> &tpl) const // tpl is a tuple of 2 elements: (count_mom, ice)
        {
          real_t count_mom = thrust::get<0>(tpl);
          int ice_flag = thrust::get<1>(tpl);
          return (ice_flag==0) ? count_mom : real_t(0);
        }
      };

      // Functor to select moment of only ice particles
      template<class real_t>
      class select_ice
      {
      public:
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, int> &tpl) const // tpl is a tuple of 2 elements: (count_mom, ice)
        {
          real_t count_mom = thrust::get<0>(tpl);
          int ice_flag = thrust::get<1>(tpl);
          return (ice_flag==1) ? count_mom : real_t(0);
        }
      };

    }


    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::cond(
      const real_t &dt,
      const real_t &RH_max,
      const bool turb_cond
    ) {   

      namespace arg = thrust::placeholders;

      thrust_device::vector<real_t> &lambda_D(tmp_device_real_cell1); // real_cell used in cond.ipp
      thrust_device::vector<real_t> &lambda_K(tmp_device_real_cell2); // real_cell used in cond.ipp

      // --- calc liquid water content before cond ---
      hskpng_sort(); 
      //thrust_device::vector<real_t> &drv(tmp_device_real_cell);
      thrust_device::vector<real_t> &drv_liq(tmp_device_real_cell);
      thrust_device::vector<real_t> &drv_ice(tmp_device_real_cell3);

      // Compute per-cell 3rd wet moment before condensation (sum of n*r^3). It is stored in count_mom
      moms_all();
      moms_calc(rw2.begin(), real_t(3./2.));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) before condensation");

      // permute-copying the result to -dm_3
      // fill with 0s if not all cells will be updated in the following transform
      if(count_n!=n_cell) {
        thrust::fill(drv_liq.begin(), drv_liq.end(), real_t(0.));
        thrust::fill(drv_ice.begin(), drv_ice.end(), real_t(0.));
      }

      // fill drv_liq with liquid moment before conde
      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(count_mom.begin(), ice.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(count_mom.begin() + count_n, ice.begin() + count_n)), // input - 1st arg
        thrust::make_permutation_iterator(drv_liq.begin(), count_ijk.begin()), // output
        detail::select_liq<real_t>()
      );

      // fill drv_ice with ice moment before cond
      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(count_mom.begin(), ice.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(count_mom.begin() + count_n, ice.begin() + count_n)), // input - 1st arg
        thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()), // output
        detail::select_ice<real_t>()
      );

      // thrust::transform(
      //   count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
      //   thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
      //   thrust::negate<real_t>()
      // );

      auto hlpr_zip_iter = thrust::make_zip_iterator(thrust::make_tuple(
        thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
        thrust::make_permutation_iterator(rv.begin(), ijk.begin()),
        thrust::make_permutation_iterator(T.begin(), ijk.begin()),
        thrust::make_permutation_iterator(eta.begin(), ijk.begin()),
        rd3.begin(),
        kpa.begin(),
        vt.begin(),
        thrust::make_permutation_iterator(lambda_D.begin(), ijk.begin()),
        thrust::make_permutation_iterator(lambda_K.begin(), ijk.begin()),
        ice.begin()
      ));

      // calculating drop growth in a timestep using backward Euler 
      // TODO: both calls almost identical, use std::bind or sth?
      if(turb_cond)
      {
        thrust_device::vector<real_t> &RH_plus_ssp(tmp_device_real_part2);
        thrust::transform(
          ssp.begin(), ssp.end(),
          thrust::make_permutation_iterator(RH.begin(), ijk.begin()),
          RH_plus_ssp.begin(),
          arg::_1 + arg::_2
        );

        // no RH_i, because we dont allow ice with turb_cond
        /*
        thrust::transform(
          rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
          thrust::make_zip_iterator(      // input - 2nd arg
            thrust::make_tuple(
              hlpr_zip_iter,
              thrust::make_permutation_iterator(p.begin(), ijk.begin()),
              RH_plus_ssp.begin(),
              thrust::make_constant_iterator<real_t>(0) // dummy RH_i just to make it compile
            )
          ), 
          rw2.begin(),                    // output
          detail::advance_rw2<real_t>(dt, RH_max)
        );
        */
      }
      else
        thrust::transform(
          rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
          thrust::make_zip_iterator(      // input - 2nd arg
            thrust::make_tuple(
              hlpr_zip_iter,
              thrust::make_permutation_iterator(p.begin(), ijk.begin()),
              thrust::make_permutation_iterator(RH.begin(), ijk.begin()),
              thrust::make_permutation_iterator(RH_i.begin(), ijk.begin())
            )
          ), 
          rw2.begin(),                    // output
          detail::advance_rw2<real_t>(dt, RH_max)
        );
      nancheck(rw2, "rw2 after condensation (no sub-steps");

      //  Compute per-cell 3rd wet moment after condensation. It is stored in count_mom
      moms_calc(rw2.begin(), real_t(3./2.));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) after condensation");

      // fill count_mom with liquid only
      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(count_mom.begin(), ice.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(count_mom.begin() + count_n, ice.begin() + count_n)), // input - 1st arg
        thrust::make_permutation_iterator(count_mom.begin(), count_ijk.begin()), // output
        detail::select_liq<real_t>()
      );

      // calculating the change in liquid moment
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv_liq.begin(), count_ijk.begin()), // input - 2nd arg
        thrust::make_permutation_iterator(drv_liq.begin(), count_ijk.begin()), // output
        thrust::minus<real_t>()
      );


      //  Compute again the per-cell 3rd wet moment after condensation. It is stored in count_mom
      moms_calc(rw2.begin(), real_t(3./2.));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) after condensation");

      // fill count_mom_ with ice only
      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(count_mom.begin(), ice.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(count_mom.begin() + count_n, ice.begin() + count_n)), // input - 1st arg
        thrust::make_permutation_iterator(count_mom.begin(), count_ijk.begin()), // output
        detail::select_ice<real_t>()
      );

      // calculating the change in ice moment
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()), // input - 2nd arg
        thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()), // output
        thrust::minus<real_t>()
      );

      // // adding the third moment after condensation to dm_3
      // thrust::transform(
      //   count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
      //   thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // input - 2nd arg
      //   thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
      //   thrust::plus<real_t>()
      // );

      // update th and rv according to changes in third specific wet moments
      update_th_rv(drv_liq);
      update_th_rv_subl(drv_ice);
    }
  };  
};
