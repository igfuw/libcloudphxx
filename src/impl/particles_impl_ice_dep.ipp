// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/iterator/transform_iterator.h>
#include <libcloudph++/common/maxwell-mason.hpp>
#include <libcloudph++/common/kappa_koehler.hpp>
#include <libcloudph++/common/kelvin_term.hpp>
#include <libcloudph++/common/transition_regime.hpp>
#include <libcloudph++/common/ventil.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {

    namespace detail
    {
      template <typename real_t>
      struct advance_rw2_minfun_ice
      {
        const quantity<si::area,              real_t> r2_old;
        const quantity<si::time,              real_t> dt;
        const quantity<si::mass_density,      real_t> rhod;
        const quantity<si::dimensionless,     real_t> rv;
        const quantity<si::temperature,       real_t> T;
        const quantity<si::pressure,          real_t> p;
        const quantity<si::dimensionless,     real_t> RH_i;
        const quantity<si::dynamic_viscosity, real_t> eta;
        const quantity<si::volume,            real_t> rd3;
        const quantity<si::dimensionless,     real_t> kpa;
        const quantity<si::velocity,          real_t> vt;
        const quantity<si::dimensionless,     real_t> RH_max;
        const quantity<si::length,            real_t> lambda_D;
        const quantity<si::length,            real_t> lambda_K;

        // ctor
        BOOST_GPU_ENABLED
        advance_rw2_minfun_ice(
          const real_t &dt,
          const real_t &rw2,
          const thrust::tuple<thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t>, real_t, real_t> &tpl,
          const real_t &RH_max
        ) :
          dt(dt * si::seconds),
          r2_old(rw2 * si::square_metres),
          rhod(    thrust::get<0>(thrust::get<0>(tpl)) * si::kilograms / si::cubic_metres),
          rv(      thrust::get<1>(thrust::get<0>(tpl))),
          T(       thrust::get<2>(thrust::get<0>(tpl)) * si::kelvins),
          eta(     thrust::get<3>(thrust::get<0>(tpl)) * si::pascals * si::seconds),
          rd3(     thrust::get<4>(thrust::get<0>(tpl)) * si::cubic_metres),
          kpa(     thrust::get<5>(thrust::get<0>(tpl))),
          vt(      thrust::get<6>(thrust::get<0>(tpl)) * si::metres_per_second),
          p(       thrust::get<1>(tpl) * si::pascals),
          RH_i(      thrust::get<2>(tpl)),
          lambda_D(thrust::get<7>(thrust::get<0>(tpl)) * si::metres),
          lambda_K(thrust::get<8>(thrust::get<0>(tpl)) * si::metres),
          RH_max(RH_max)
        {}

        BOOST_GPU_ENABLED
        quantity<divide_typeof_helper<si::area, si::time>::type, real_t> drw2_dt(const quantity<si::area, real_t> &rw2) const
        {
          using namespace common::maxwell_mason;
          using namespace common::kappa_koehler;
          using namespace common::kelvin;
          using common::moist_air::D_0;
          using common::moist_air::K_0;
          using common::moist_air::c_pd;
          using common::transition_regime::beta;
          using common::ventil::Sh;
          using common::ventil::Nu;
#if !defined(__NVCC__)
          using std::sqrt;
#endif

          const quantity<si::length, real_t> rw  = sqrt(real_t(rw2 / si::square_metres)) * si::metres;
          const quantity<si::volume, real_t> rw3 = rw * rw * rw;;

          const quantity<si::dimensionless, real_t>
            Re = common::ventil::Re(vt, rw, rhod, eta),
            Sc = common::ventil::Sc(eta, rhod, D_0<real_t>()), // TODO? cache
            Pr = common::ventil::Pr(eta, c_pd<real_t>(), K_0<real_t>()); // TODO? cache

          const quantity<common::diffusivity, real_t>
            D = D_0<real_t>() * beta(lambda_D / rw) * (Sh(Sc, Re) / 2);

          const quantity<common::thermal_conductivity, real_t>
            K = K_0<real_t>() * beta(lambda_K / rw) * (Nu(Pr, Re) / 2);

          return real_t(2) * rdrdt_i(
            D,
            K,
            rhod * rv,
            T,
            p,
            RH_i > RH_max ? RH_max : RH_i
          );
        }

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2_unitless) const
        {
          const quantity<si::area, real_t> rw2 = rw2_unitless * si::square_metres;
          return (r2_old + dt * drw2_dt(rw2) - rw2) / si::square_metres;
        }
      };


      template <typename real_t>
      struct advance_ice_ac
      {
        const real_t dt, RH_max;

        advance_ice_ac(const real_t &dt, const real_t &RH_max)
          : dt(dt), RH_max(RH_max) {}

        BOOST_GPU_ENABLED
        thrust::tuple<real_t, real_t> operator()(
            const thrust::tuple<real_t, real_t> &ac_old,
            const thrust::tuple<
                thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t>,
                real_t, real_t> &tpl
        ) const
        {
#if !defined(__NVCC__)
          using std::max;
          using std::isnan;
          using std::isinf;
#endif
          const real_t a_old = thrust::get<0>(ac_old);
          const real_t c_old = thrust::get<1>(ac_old);

          // Skip liquid droplets
          if (a_old <= 0 || c_old <= 0)
            return ac_old;

          advance_rw2_minfun_ice<real_t> f_a(dt, a_old * a_old, tpl, RH_max);
          advance_rw2_minfun_ice<real_t> f_c(dt, c_old * c_old, tpl, RH_max);

          const real_t da_dt = (f_a.drw2_dt(a_old * a_old * si::square_metres) / (2 * a_old * si::metres))
                                * si::seconds / si::metres;
          const real_t dc_dt = (f_c.drw2_dt(c_old * c_old * si::square_metres) / (2 * c_old * si::metres))
                                * si::seconds / si::metres;

          // forward Euler for simplicity
          const real_t a_new = max(a_old + dt * da_dt, real_t(1e-9));
          const real_t c_new = max(c_old + dt * dc_dt, real_t(1e-9));

          return thrust::make_tuple(a_new, c_new);
        }
      };

    }


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
      auto drv_ice_g = tmp_device_real_cell.get_guard();
      thrust_device::vector<real_t> &drv_ice = drv_ice_g.get();
      if(step == 0)
        reset_guardp(ice_mass_gp, tmp_device_real_cell);
      thrust_device::vector<real_t> &ice_mass = ice_mass_gp->get();

      // Compute per-cell 3rd moment of ice before deposition. It is stored in count_mom
      if(step == 0)
      {
        moms_gt0(ice_a.begin()); // choose ice particles (ice_a>0)
        moms_calc(thrust::make_transform_iterator(
          thrust::make_zip_iterator(thrust::make_tuple(ice_a.begin(), ice_c.begin(), ice_rho.begin())),
          detail::ice_mass<real_t>()
          ),
          real_t(1));
        nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd ice moment) before deposition");

        // fill with 0s if not all cells have particles
        if(count_n!=n_cell) {
          thrust::fill(drv_ice.begin(), drv_ice.end(), real_t(0.));
          thrust::fill(ice_mass.begin(), ice_mass.end(), real_t(0.));
        }

        thrust::transform(
          count_mom.begin(), count_mom.begin() + count_n,
          thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()),
          thrust::negate<real_t>()
        );
      }
      else // copy ice_mass from previous step
      {
        // drv = -ice_mass precond
        thrust::transform(
          ice_mass.begin(), ice_mass.end(),
          drv_ice.begin(),
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
        detail::advance_ice_ac<real_t>(dt, RH_max)
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

      // Adding the ice volume after deposition to drv_ice
      if(step < sstp_cond - 1)
      {
        thrust::copy(
          count_mom.begin(), count_mom.begin() + count_n,                        // input - 1st arg
          thrust::make_permutation_iterator(ice_mass.begin(), count_ijk.begin())  // output
        );

        // adding the third moment after deposition to ice_mass
        thrust::transform(
          ice_mass.begin(), ice_mass.end(),
          drv_ice.begin(),
          drv_ice.begin(),
          thrust::plus<real_t>()
        );
      }
      else // last step, calculate change in 3rd moment and update th and rv
      {
        thrust::transform(
          count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
          thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()), // input - 2nd arg
          thrust::make_permutation_iterator(drv_ice.begin(), count_ijk.begin()), // output
          thrust::plus<real_t>()
        );
        ice_mass_gp.reset(); // destroy guard to tmp array that stored ice_mass
      }

      // update th and rv according to change in third specific moment
      update_th_rv(drv_ice, impl::phase_change::deposition);
      }
    };
  }