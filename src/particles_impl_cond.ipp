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
#include <libcloudph++/common/detail/bisect.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename n_t, typename real_t>
      struct dm_3_summator
      {   
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<n_t, real_t> &tpl) const
        {
          const n_t n = thrust::get<0>(tpl);
          const real_t rw2 = thrust::get<1>(tpl);
          return n * pow(rw2, real_t(3./2));
        }
      };  
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::cond_dm3_helper() // TODO: move it into a common_count_mom?
    {   
      // TODO... or at least these typedefs!
      typedef thrust::permutation_iterator<
        typename thrust_device::vector<n_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_n_t;
      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_r_t;
      typedef thrust::zip_iterator<thrust::tuple<pi_n_t, pi_r_t> > zip_it_t;

      thrust::pair<
	typename thrust_device::vector<thrust_size_t>::iterator,
	typename thrust_device::vector<real_t>::iterator
      > n = thrust::reduce_by_key(
	sorted_ijk.begin(), sorted_ijk.end(),   // input - keys
	thrust::transform_iterator<             // input - values
          detail::dm_3_summator<n_t, real_t>,
          zip_it_t,
          real_t
        >( 
	  zip_it_t(thrust::make_tuple(
	    pi_n_t(this->n.begin(), sorted_id.begin()),
	    pi_r_t(rw2.begin(),     sorted_id.begin())
	  )),
	  detail::dm_3_summator<n_t, real_t>()
	),
	count_ijk.begin(),                      // output - keys
	count_mom.begin()                       // output - values
      );  
      count_n = n.first - count_ijk.begin();
      assert(count_n > 0 && count_n <= n_cell);
    }


    namespace detail
    {
      template <typename real_t>
      struct dth
      {
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::dimensionless, real_t> 
            drv      = thrust::get<0>(tpl);
          const quantity<si::temperature, real_t> 
            T        = thrust::get<1>(tpl) * si::kelvins;
          const quantity<si::temperature, real_t> 
            th       = thrust::get<2>(tpl) * si::kelvins;

          return drv * common::theta_dry::d_th_d_rv(T, th) / si::kelvins;
        }
      };

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

          // ''predictor'' step using explicit Euler scheme
          const advance_rw2_minfun<real_t> f(dt, rw2_old, tpl, RH_max); 
          const real_t drw2 = dt * f.drw2_dt(rw2_old * si::square_metres) * si::seconds / si::square_metres;

          if (drw2 == 0) return rw2_old;

          // ''corrector'' step using implicit Euler scheme
          const real_t rd2 = pow(thrust::get<6>(tpl), real_t(2./3));
          const real_t tol_r2 = rd2 / 10; // TODO !!! (think of a better value, document)

          const real_t mlt = 2; // results in explicit Euler scheme if root not found // TODO: investidate why not found...
          
          const real_t 
            a = max(rd2, rw2_old + min(real_t(0), mlt * drw2)),
            b =          rw2_old + max(real_t(0), mlt * drw2);
 
          // numerics (drw2 != 0 but a==b)
          if (a == b) return rw2_old;

          if (drw2 > 0) 
            return common::detail::bisect(f, a, b, tol_r2, drw2); // for implicit Euler its equal to min_fun(x_old) 
          else
            return common::detail::bisect(f, a, b, tol_r2, f(a), drw2);
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

      //
      thrust_device::vector<real_t> &drv(tmp_device_real_cell);

      // calculating the 3rd wet moment before condensation (still not divided by dv)
      cond_dm3_helper();

      // permute-copying the result to -dm_3
      thrust::fill(drv.begin(), drv.end(), 0);
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
	    thrust::make_permutation_iterator(rhod.begin(),    ijk.begin()),
	    thrust::make_permutation_iterator(rv.begin(),      ijk.begin()),
	    thrust::make_permutation_iterator(T.begin(),       ijk.begin()),
	    thrust::make_permutation_iterator(p.begin(),       ijk.begin()),
	    thrust::make_permutation_iterator(RH.begin(),      ijk.begin()),
	    thrust::make_permutation_iterator(eta.begin(),     ijk.begin()),
	    rd3.begin(),
	    kpa.begin(),
            vt.begin()
          )
        ), 
	rw2.begin(),                    // output
        detail::advance_rw2<real_t>(dt, RH_max)
      );

      // calculating the 3rd wet moment after condensation (still not divided by dv)
      cond_dm3_helper();

      // adding the third moment after condensation to dm_3
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // input - 2nd arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );

      // multiplying dv*dm_3 by -rho_w*4/3*pi
      thrust::transform(
        drv.begin(), drv.end(),                  // input - 1st arg
        thrust::make_constant_iterator<real_t>(  // input - 2nd arg
          - common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres
          * real_t(4./3) * pi<real_t>()
        ),
        drv.begin(),                             // output
        thrust::multiplies<real_t>()
      );

      // dividing by dv
      thrust::transform(
        drv.begin(), drv.end(),  // input - 1st arg
        dv.begin(),              // input - 2nd arg
        drv.begin(),             // output
        thrust::divides<real_t>() 
      ); 

      // dividing d(rhod_rv) by rhod
      thrust::transform(
        drv.begin(), drv.end(),  // input - 1st arg
        rhod.begin(),            // input - 2nd arg
        drv.begin(),             // output (in place)
        thrust::divides<real_t>()
      );

      // updating rv 
      assert(*thrust::min_element(rv.begin(), rv.end()) >= 0);
      thrust::transform(
        rv.begin(), rv.end(),  // input - 1st arg
        drv.begin(),           // input - 2nd arg
        rv.begin(),            // output
        thrust::plus<real_t>() 
      );
      assert(*thrust::min_element(rv.begin(), rv.end()) >= 0);

      // updating th
      {
        typedef thrust::zip_iterator<thrust::tuple<
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator
        > > zip_it_t;
 
	thrust::transform(
	  th.begin(), th.end(),          // input - 1st arg
	  thrust::transform_iterator<    // input - 2nd arg
	    detail::dth<real_t>,
	    zip_it_t,
	    real_t
	  >(
            zip_it_t(thrust::make_tuple(  // args (note: rv cannot be used here as already modified)
	      drv.begin(),      // 
	      T.begin(),        // dth = drv * d_th_d_rv(T, th)
	      th.begin()        //
	    )),
	    detail::dth<real_t>()     // func
	  ),
	  th.begin(),                 // output
	  thrust::plus<real_t>()
	);
      }
    }
  };  
};
