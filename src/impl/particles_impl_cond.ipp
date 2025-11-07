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
      const real_t &RH_max,
      const bool turb_cond,
      const int step
    ) {   

      namespace arg = thrust::placeholders;

      thrust_device::vector<real_t> &lambda_D(lambda_D_gp->get());
      thrust_device::vector<real_t> &lambda_K(lambda_K_gp->get());

      hskpng_sort();

      // Vector to store 3rd moment
      auto drv_liq_g = tmp_device_real_cell.get_guard();
      thrust_device::vector<real_t> &drv_liq = drv_liq_g.get();
      if(step == 0)
        reset_guardp(rw_mom3_gp, tmp_device_real_cell);
      thrust_device::vector<real_t> &rw_mom3 = rw_mom3_gp->get();

      // Compute per-cell 3rd moment of liquid droplets before condensation. It is stored in count_mom
      if(step == 0)
      {
        moms_gt0(rw2.begin()); // choose liquid particles (rw2>0)
        moms_calc(rw2.begin(), real_t(3./2.));
        nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) before condensation");

        // fill with 0s if not all cells have particles
        if(count_n!=n_cell) {
          thrust::fill(drv_liq.begin(), drv_liq.end(), real_t(0.));
          thrust::fill(rw_mom3.begin(), rw_mom3.end(), real_t(0.));
        }

        thrust::transform(
          count_mom.begin(), count_mom.begin() + count_n,
          thrust::make_permutation_iterator(drv_liq.begin(), count_ijk.begin()),
          thrust::negate<real_t>()
        );
      }
      else // copy rw_mom3 from previous step
      {
        // drv = -rw_mom3 precond
        thrust::transform(
          rw_mom3.begin(), rw_mom3.end(),
          drv_liq.begin(),
          thrust::negate<real_t>()
        );
      }

      auto hlpr_zip_iter = thrust::make_zip_iterator(thrust::make_tuple(
        thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
        thrust::make_permutation_iterator(rv.begin(), ijk.begin()),
        thrust::make_permutation_iterator(T.begin(), ijk.begin()),
        thrust::make_permutation_iterator(eta.begin(), ijk.begin()),
        rd3.begin(),
        kpa.begin(),
        vt.begin(),
        thrust::make_permutation_iterator(lambda_D.begin(), ijk.begin()),
        thrust::make_permutation_iterator(lambda_K.begin(), ijk.begin())
      ));

      // calculating drop growth in a timestep using backward Euler 
      // TODO: both calls almost identical, use std::bind or sth?
      if(turb_cond)
      {
        auto RH_plus_ssp_g = tmp_device_real_part.get_guard();
        thrust_device::vector<real_t> &RH_plus_ssp = RH_plus_ssp_g.get();
        thrust::transform(
          ssp.begin(), ssp.end(),
          thrust::make_permutation_iterator(RH.begin(), ijk.begin()),
          RH_plus_ssp.begin(),
          arg::_1 + arg::_2
        );

        // condensation for liquid droplets
        thrust::transform(
          rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
          thrust::make_zip_iterator(      // input - 2nd arg
            thrust::make_tuple(
              hlpr_zip_iter,
              thrust::make_permutation_iterator(p.begin(), ijk.begin()),
              RH_plus_ssp.begin()
            )
          ), 
          rw2.begin(),                    // output
          detail::advance_rw2<real_t>(dt, RH_max)
        );
      }
      else
      {
        // condensation for liquid droplets
        thrust::transform(
          rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
          thrust::make_zip_iterator(      // input - 2nd arg
            thrust::make_tuple(
              hlpr_zip_iter,
              thrust::make_permutation_iterator(p.begin(), ijk.begin()),
              thrust::make_permutation_iterator(RH.begin(), ijk.begin())
            )
          ), 
          rw2.begin(),                    // output
          detail::advance_rw2<real_t>(dt, RH_max)
        );
        nancheck(rw2, "rw2 after condensation (no sub-steps");

        // Compute per-cell 3rd moment of liquid droplets after condensation. It is stored in count_mom
        moms_gt0(rw2.begin()); // choose liquid particles (rw2>0)
        moms_calc(rw2.begin(), real_t(3./2.));
        nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) after condensation");

        // Adding the third liquid moment after condensation to drv_liq
        if(step < sstp_cond - 1)
        {
          thrust::copy(
            count_mom.begin(), count_mom.begin() + count_n,                        // input - 1st arg
            thrust::make_permutation_iterator(rw_mom3.begin(), count_ijk.begin())  // output
          );

          // adding the third moment after condensation to dm_3
          thrust::transform(
            rw_mom3.begin(), rw_mom3.end(),
            drv_liq.begin(),
            drv_liq.begin(),
            thrust::plus<real_t>()
          );
        }
        else // last step, calculate change in 3rd moment and update th and rv
        {
          thrust::transform(
            count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
            thrust::make_permutation_iterator(drv_liq.begin(), count_ijk.begin()), // input - 2nd arg
            thrust::make_permutation_iterator(drv_liq.begin(), count_ijk.begin()), // output
            thrust::plus<real_t>()
          );
          rw_mom3_gp.reset(); // destroy guard to tmp array that stored 3rd moment of rw
        }

        // update th and rv according to change in third specific wet moment
        update_th_rv(drv_liq, impl::phase_change::condensation);


        if (opts_init.ice_switch)
        {
          // Vector to store 3rd moment
          auto drv_ice_g = tmp_device_real_cell.get_guard();
          thrust_device::vector<real_t> &drv_ice = drv_ice_g.get();
          if(step == 0)
            reset_guardp(rw_mom3_gp, tmp_device_real_cell);
          thrust_device::vector<real_t> &rw_mom3 = rw_mom3_gp->get();

          // Compute per-cell 3rd moment of ice before deposition. It is stored in count_mom
          if(step == 0)
          {
            moms_gt0(ice_a.begin()); // choose ice particles (ice_a>0)
            moms_calc(thrust::make_transform_iterator(
              thrust::make_zip_iterator(thrust::make_tuple(ice_a.begin(), ice_c.begin())), detail::ice_vol<real_t>()
              ),
              real_t(1));
            nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd ice moment) before deposition");
            if(count_n!=n_cell) {
              thrust::fill(drv_ice.begin(), drv_ice.end(), real_t(0.));
            }

            thrust::transform(
              count_mom.begin(), count_mom.begin() + count_n,
              thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()),
              thrust::negate<real_t>()
            );
          }
          else // copy rw_mom3 from previous step
          {
            // drv = -rw_mom3 precond
            thrust::transform(
              rw_mom3.begin(), rw_mom3.end(),
              drv_ice.begin(),
              thrust::negate<real_t>()
            );
          }

          // deposition for ice crystals
          thrust::transform(
            thrust::make_zip_iterator(
              thrust::make_tuple(
                ice_a.begin(),
                ice_c.begin()
              )
            ),
            thrust::make_zip_iterator(
            thrust::make_tuple(
              ice_a.end(),
              ice_c.end()
            )
          ),
            thrust::make_zip_iterator(
              thrust::make_tuple(
                hlpr_zip_iter,
                thrust::make_permutation_iterator(p.begin(), ijk.begin()),
                thrust::make_permutation_iterator(RH_i.begin(), ijk.begin())
              )
            ),
            thrust::make_zip_iterator(
              thrust::make_tuple(
                ice_a.begin(),
                ice_c.begin()
              )
            ),
            detail::advance_ice_ac<real_t>(dt, RH_max)
          );
          nancheck(ice_a, "ice_a after deposition (no sub-steps");
          nancheck(ice_c, "ice_c after deposition (no sub-steps");

          // Compute per-cell 3rd moment of ice after deposition. It is stored in count_mom
          moms_gt0(ice_a.begin()); // choose ice particles (ice_a>0)
          moms_calc(thrust::make_transform_iterator(
            thrust::make_zip_iterator(thrust::make_tuple(ice_a.begin(), ice_c.begin())), detail::ice_vol<real_t>()
            ),
            real_t(1));
          nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd ice moment) after deposition");

          // Adding the third ice moment after deposition to drv_ice
          thrust::transform(
            count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
            thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()), // input - 2nd arg
            thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()), // output
            thrust::plus<real_t>()
          );

          // update th and rv according to change in third specific wet moment
          update_th_rv(drv_ice, impl::phase_change::sublimation);
        }
      }
    };
  }
}