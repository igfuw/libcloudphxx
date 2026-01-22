// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

// #include <thrust/iterator/transform_iterator.h>
// #include <libcloudph++/common/maxwell-mason.hpp>
// #include <libcloudph++/common/kappa_koehler.hpp>
// #include <libcloudph++/common/kelvin_term.hpp>
// #include <libcloudph++/common/transition_regime.hpp>
// #include <libcloudph++/common/ventil.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::ice_dep(
      const real_t &dt,
      const real_t &RH_max,
      const int step
    ) {   

      namespace arg = thrust::placeholders;

      thrust_device::vector<real_t> &lambda_D(lambda_D_gp->get());
      thrust_device::vector<real_t> &lambda_K(lambda_K_gp->get());

      hskpng_sort();

      // Vector to store 3rd moment
      // auto d_ice_mass_g = tmp_device_real_cell.get_guard();
      // thrust_device::vector<real_t> &d_ice_mass = d_ice_mass_g.get();
      if(step == 0)
        reset_guardp(ice_mass_gp, tmp_device_real_cell);
      thrust_device::vector<real_t> &ice_mass = ice_mass_gp->get();

      if(step == 0)
      {
        if(count_n!=n_cell)  // NOTE: we rely that moms_calc was called recently (e.g. in save_liq_ice_content_before_change)
          thrust::fill(ice_mass.begin(), ice_mass.end(), real_t(0.));        
      }
      else // copy ice_mass from previous step
      {
        reset_guardp(d_ice_mass_gp, tmp_device_real_cell);
        thrust_device::vector<real_t> &d_ice_mass = d_ice_mass_gp->get();
        // drv = -ice_mass precond
        thrust::transform(
          ice_mass.begin(), ice_mass.end(),
          d_ice_mass.begin(),
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
        detail::advance_ice_ac<real_t>(dt / sstp_cond, RH_max)
      );
      nancheck(ice_a, "ice_a after deposition (no sub-steps");
      nancheck(ice_c, "ice_c after deposition (no sub-steps");

      // Compute per-cell 3rd moment of ice after deposition. It is stored in count_mom
      moms_gt0(ice_a.begin()); // choose ice particles (ice_a>0)
      moms_calc(thrust::make_transform_iterator(
        thrust::make_zip_iterator(thrust::make_tuple(ice_a.begin(), ice_c.begin(), ice_rho.begin())),
        detail::ice_mass<real_t>()
        ),
        real_t(1));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd ice moment) after deposition");

      // Adding the ice volume after deposition to d_ice_mass
      if(step < sstp_cond - 1)
      {
        thrust_device::vector<real_t> &d_ice_mass = d_ice_mass_gp->get();
        thrust::copy(
          count_mom.begin(), count_mom.begin() + count_n,                        // input - 1st arg
          thrust::make_permutation_iterator(ice_mass.begin(), count_ijk.begin())  // output
        );

        // adding the third moment after deposition to ice_mass
        thrust::transform(
          ice_mass.begin(), ice_mass.end(),
          d_ice_mass.begin(),
          d_ice_mass.begin(),
          thrust::plus<real_t>()
        );
      }
      else // last step, calculate change in 3rd moment and update th and rv
      {
        thrust_device::vector<real_t> &d_ice_mass = d_ice_mass_gp->get();
        
        thrust::transform(
          count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
          thrust::make_permutation_iterator(d_ice_mass.begin(), count_ijk.begin()), // input - 2nd arg
          thrust::make_permutation_iterator(d_ice_mass.begin(), count_ijk.begin()), // output
          thrust::plus<real_t>()
        );
        ice_mass_gp.reset(); // destroy guard to tmp array that stored ice_mass
      }

      // update th and rv according to change in third specific moment
      // update_th_rv(d_ice_mass, impl::phase_change::deposition);
      }
    };
  }