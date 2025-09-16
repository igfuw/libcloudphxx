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

      template<class real_t>
      class ice_vol
      {
      public:
        using result_type = real_t;
        template <typename Tuple>
        BOOST_GPU_ENABLED
        real_t operator()(Tuple const &tpl) const
        {
          real_t a = thrust::get<0>(tpl);
          real_t c = thrust::get<1>(tpl);
          return a * a * c;
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

      hskpng_sort();

      // Vectors to store 3rd moments
      thrust_device::vector<real_t> &drv_liq(tmp_device_real_cell);
      thrust_device::vector<real_t> &drv_ice(tmp_device_real_cell3);

      // Compute per-cell 3rd moment of liquid droplets before condensation. It is stored in count_mom
      moms_eq0(ice.begin()); // choose particles with ice=0
      moms_calc(rw2.begin(), real_t(3./2.));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) before condensation");
      if(count_n!=n_cell) {
        thrust::fill(drv_liq.begin(), drv_liq.end(), real_t(0.));
      }

      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,
        thrust::make_permutation_iterator(drv_liq.begin(), count_ijk.begin()),
        thrust::negate<real_t>()
      );

      // Compute per-cell 3rd moment of ice before sublimation. It is stored in count_mom
      moms_gt0(ice.begin()); // choose particles with ice=1
      moms_calc(thrust::make_transform_iterator(
        thrust::make_zip_iterator(thrust::make_tuple(ice_a.begin(), ice_c.begin())), detail::ice_vol<real_t>()
        ),
        real_t(1));

      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd ice moment) before sublimation");
      if(count_n!=n_cell) {
        thrust::fill(drv_ice.begin(), drv_ice.end(), real_t(0.));
      }

      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,
        thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()),
        thrust::negate<real_t>()
      );


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

      // Compute per-cell 3rd moment of liquid droplets after condensation. It is stored in count_mom
      moms_eq0(ice.begin()); // choose particles with ice=0
      moms_calc(rw2.begin(), real_t(3./2.));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) after condensation");

      // Adding the third liquid moment after condensation to drv_liq
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv_liq.begin(), count_ijk.begin()), // input - 2nd arg
        thrust::make_permutation_iterator(drv_liq.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );

      // Compute per-cell 3rd moment of ice after sublimation. It is stored in count_mom
      moms_gt0(ice.begin()); // choose particles with ice=1
      moms_calc(thrust::make_transform_iterator(
        thrust::make_zip_iterator(thrust::make_tuple(ice_a.begin(), ice_c.begin())), detail::ice_vol<real_t>()
        ),
        real_t(1));

      // Adding the third ice moment after sublimation to drv_ice
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()), // input - 2nd arg
        thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );

      // update th and rv according to changes in third specific wet moments
      update_th_rv(drv_liq, impl::phase_change::condensation);
      update_th_rv(drv_ice, impl::phase_change::sublimation);
    }
  };
};