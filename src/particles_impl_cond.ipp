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
#include <libcloudph++/common/mean_free_path.hpp>
#include <libcloudph++/common/detail/toms748.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct advance_rw2_minfun
      {
        const quantity<si::area,              real_t> rw2_old;
        const quantity<si::time,              real_t> dt;
	const quantity<si::mass_density,      real_t> rhod;
	const quantity<si::dimensionless,     real_t> rv;
	const quantity<si::temperature,       real_t> T;
	const quantity<si::pressure,          real_t> p;
	const quantity<si::dimensionless,     real_t> RH;
	const quantity<si::dynamic_viscosity, real_t> eta;
	const quantity<si::volume,            real_t> rd3;
	const quantity<si::dimensionless,     real_t> kpa;
	const quantity<si::velocity,          real_t> vt;
        const quantity<si::dimensionless,     real_t> RH_max;

        // ctor
        BOOST_GPU_ENABLED
        advance_rw2_minfun(
          const real_t &dt,
          const real_t &rw2,
          const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t> &tpl,
          const real_t &RH_max
        ) : 
          dt(dt * si::seconds), 
          rw2_old(rw2 * si::square_metres),
          rhod(    thrust::get<0>(tpl) * si::kilograms / si::cubic_metres),
          rv(      thrust::get<1>(tpl)),
          T(       thrust::get<2>(tpl) * si::kelvins),
          p(       thrust::get<3>(tpl) * si::pascals),
          RH(      thrust::get<4>(tpl)),
          eta(     thrust::get<5>(tpl) * si::pascals * si::seconds),
          rd3(     thrust::get<6>(tpl) * si::cubic_metres),
          kpa(     thrust::get<7>(tpl)),
          vt(      thrust::get<8>(tpl) * si::metres_per_second),
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
          using common::mean_free_path::lambda_D;
          using common::mean_free_path::lambda_K;
          using common::ventil::Sh;
          using common::ventil::Nu;
          using std::sqrt;

	  const quantity<si::length, real_t> rw  = sqrt(real_t(rw2 / si::square_metres)) * si::metres; 
	  const quantity<si::volume, real_t> rw3 = rw * rw * rw;;

          // TODO: common::moist_air:: below should not be needed
          // TODO: ventilation as option
          const quantity<si::dimensionless, real_t>
            Re = common::ventil::Re(vt, rw, rhod, eta), 
            Sc = common::ventil::Sc(eta, rhod, D_0<real_t>()), // TODO? cache
            Pr = common::ventil::Pr(eta, c_pd<real_t>(), K_0<real_t>()); // TODO? cache

          const quantity<common::diffusivity, real_t> 
            D = D_0<real_t>() * beta(lambda_D(T)    / rw) * (Sh(Sc, Re) / 2); // TODO: cache lambdas

          const quantity<common::thermal_conductivity, real_t> 
            K = K_0<real_t>() * beta(lambda_K(T, p) / rw) * (Nu(Pr, Re) / 2);

          return real_t(2) * rdrdt( 
            D,
            K,
            rhod * rv, 
            T, 
            p, 
            RH > RH_max ? RH_max : RH, 
            a_w(rw3, rd3, kpa),
            klvntrm(rw, T)
          );
        }

        // backward Euler scheme:
	// rw2_new = rw2_old + f_rw2(rw2_new) * dt
	// rw2_new = rw2_old + 2 * rw * f_rw(rw2_new) * dt
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2_unitless) const
        {
	  const quantity<si::area, real_t> rw2 = rw2_unitless * si::square_metres; 
          return (rw2_old + dt * drw2_dt(rw2) - rw2) / si::square_metres;
        }
      };

      template <typename real_t>
      struct advance_rw2
      {
        const real_t dt, RH_max;

        advance_rw2(const real_t &dt, const real_t &RH_max) : dt(dt), RH_max(RH_max) {}

        BOOST_GPU_ENABLED
        real_t operator()(
          const real_t &rw2_old, 
          const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t> &tpl
        ) const {
#if !defined(__NVCC__)
	  using std::min;
	  using std::max;
	  using std::pow;
	  using std::abs;
#endif

          const advance_rw2_minfun<real_t> f(dt, rw2_old, tpl, RH_max); 
          const real_t drw2 = dt * f.drw2_dt(rw2_old * si::square_metres) * si::seconds / si::square_metres;

          if (drw2 == 0) return rw2_old;

          const real_t rd2 = pow(thrust::get<6>(tpl), real_t(2./3));
          
          const int mlt = 2; // arbitrary!
 
          const real_t 
            a = max(rd2, rw2_old + min(real_t(0), mlt * drw2)),
            b =          rw2_old + max(real_t(0), mlt * drw2);

          // numerics (drw2 != 0 but a==b)
          if (a == b) return rw2_old;

          real_t fa, fb;

          if (drw2 > 0) 
          {
            fa = drw2; // for implicit Euler its equal to min_fun(x_old) 
            fb = f(b);
          }
          else
          {
            fa = f(a);
            fb = drw2; // for implicit Euler its equal to min_fun(x_old) 
          }

          // root-finding ill posed => explicit Euler 
	  if (fa * fb > 0) return rw2_old + drw2;

          // otherwise implicit Euler
	  return common::detail::toms748_solve(f, a, b, fa, fb); 
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::cond(
      const real_t &dt,
      const real_t &RH_max
    ) {   
      // prerequisite
      hskpng_sort(); 
      thrust_device::vector<real_t> &drv(tmp_device_real_cell);

      // calculating the 3rd wet moment before condensation (still not divided by dv)
      moms_calc_cond(rw2.begin(), real_t(3./2.));

      // permute-copying the result to -dm_3
      // fill with 0s if not all cells will be updated in the following transform
      if(count_n!=n_cell)  thrust::fill(drv.begin(), drv.end(), real_t(0.));
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::negate<real_t>()
      );

      // same vectors as in sstp_step
      thrust_device::vector<real_t> &Tp(tmp_device_real_part3),
                                     pp(tmp_device_real_part4),
                                     RHp(tmp_device_real_part5),
                                     etap(tmp_device_real_part6);

      // calculating drop growth in a timestep using backward Euler 
      thrust::transform(
        rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
        thrust::make_zip_iterator(      // input - 2nd arg
          thrust::make_tuple(
      sstp_tmp_rh.begin(),
      sstp_tmp_rv.begin(),
      Tp.begin(),
      pp.begin(),
      RHp.begin(),
      etap.begin(),
	    rd3.begin(),
	    kpa.begin(),
            vt.begin()
          )
        ), 
	rw2.begin(),                    // output
        detail::advance_rw2<real_t>(dt, RH_max)
      );

      // calculating the 3rd wet moment after condensation (still not divided by dv)
      moms_calc_cond(rw2.begin(), real_t(3./2.));

      // adding the third moment after condensation to dm_3
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // input - 2nd arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );

      // update th and rv according to changes in third wet moment
      update_th_rv(drv, true);
    }
  };  
};
