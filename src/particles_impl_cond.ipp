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

    template <typename real_t, int device>
    void particles<real_t, device>::impl::cond_dm3_helper() // TODO: move it into a common_count_mom?
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
      struct drhod_th
      {
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t> &tpl) const
        {
          const quantity<divide_typeof_helper<si::mass_density, si::time >::type, real_t> 
            drhod_rv = thrust::get<0>(tpl) * si::kilograms / si::cubic_metres / si::seconds;
          const quantity<si::mass_density, real_t> 
            rhod     = thrust::get<1>(tpl) * si::kilograms / si::cubic_metres;
          const quantity<si::temperature, real_t> 
            T        = thrust::get<2>(tpl) * si::kelvins;
          const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> 
            rhod_th  = thrust::get<3>(tpl) * si::kilograms / si::cubic_metres * si::kelvins;

          return 
            drhod_rv 
            / rhod 
            * common::theta_dry::d_rhodtheta_d_rv(T, rhod_th) 
            / si::kilograms * si::cubic_metres / si::kelvins * si::seconds;
        }
      };

      template <typename real_t>
      struct advance_rw2_minfun
      {
        const quantity<si::area,          real_t> rw2_old;
        const quantity<si::time,          real_t> dt;
	const quantity<si::mass_density,  real_t> rhod_rv;
	const quantity<si::temperature,   real_t> T;
	const quantity<si::pressure,      real_t> p;
	const quantity<si::dimensionless, real_t> RH;
	const quantity<si::volume,        real_t> rd3;
	const quantity<si::dimensionless, real_t> kpa;

const quantity<si::dimensionless, real_t> RH_max = 1.01; // TODO!

        // ctor
        advance_rw2_minfun(
          const real_t &dt,
          const real_t &rw2,
          const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t> &tpl
        ) : 
          dt(dt * si::seconds), 
          rw2_old(rw2 * si::square_metres),
          rhod_rv(thrust::get<0>(tpl) * si::kilograms / si::cubic_metres),
          T(thrust::get<1>(tpl) * si::kelvins),
          p(thrust::get<2>(tpl) * si::pascals),
          RH(thrust::get<3>(tpl)),
          rd3(thrust::get<4>(tpl) * si::cubic_metres),
          kpa(thrust::get<5>(tpl))
        {}

        BOOST_GPU_ENABLED
        quantity<divide_typeof_helper<si::area, si::time>::type, real_t> drw2_dt(const quantity<si::area, real_t> &rw2) const
        {
          using namespace common::maxwell_mason;
          using namespace common::kappa_koehler;
          using namespace common::kelvin;

using std::sqrt;

	  const quantity<si::length, real_t> rw  = sqrt(real_t(rw2 / si::square_metres)) * si::metres; 
	  const quantity<si::volume, real_t> rw3 = rw * rw * rw;;

          return real_t(2) * rdrdt( 
            rhod_rv, T, p, RH > RH_max ? RH_max : RH, 
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
        const real_t dt;
        advance_rw2(const real_t &dt) : dt(dt) {}
        real_t operator()(
          const real_t &rw2_old, 
          const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t> &tpl
        ) const {
#if !defined(__NVCC__)
	  using std::min;
	  using std::max;
	  using std::pow;
	  using std::abs;
#endif

          // ''predictor'' step using explicit Euler scheme
          const advance_rw2_minfun<real_t> f(dt, rw2_old, tpl);
          const real_t drw2 = dt * f.drw2_dt(rw2_old * si::square_metres) * si::seconds / si::square_metres;

          if (drw2 == 0) return rw2_old;

          // ''corrector'' step using implicit Euler scheme
          const real_t rd2 = pow(thrust::get<4>(tpl), real_t(2./3));
          const real_t tol_r2 = rd2 / 10; // TODO !!!

          const real_t mlt = 2; // results in explicit Euler scheme if root not found // TODO: investidate why not found...
          
          const real_t 
            a = max(rd2, rw2_old + min(real_t(0), mlt * drw2)),
            b =          rw2_old + max(real_t(0), mlt * drw2);
 
          if (drw2 > 0) 
            return common::detail::bisect(f, a, b, tol_r2, drw2); // for implicit Euler its equal to min_fun(x_old) 
          else
            return common::detail::bisect(f, a, b, tol_r2, f(a), drw2);
        }
      };
    };

    template <typename real_t, int device>
    void particles<real_t, device>::impl::cond(real_t dt)
    {   
      // prerequisites
      hskpng_sort(); // TODO: the same with T,p,r,RH? (and dependencies among T,p,r,RH!)

      //
      thrust_device::vector<real_t> &drhod_rv(tmp_device_real_cell);

      // calculating the 3rd wet moment before condensation (still not divided by dv)
      cond_dm3_helper();

      // permute-copying the result to -dm_3
      thrust::fill(drhod_rv.begin(), drhod_rv.end(), 0);
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,            // input - 1st arg
        thrust::make_permutation_iterator(drhod_rv.begin(), count_ijk.begin()), // output
        thrust::negate<real_t>()
      );

      // calculating drop growth in a timestep using backward Euler 
      thrust::transform(
        rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
        thrust::make_zip_iterator(      // input - 2nd arg
          thrust::make_tuple(
	    thrust::make_permutation_iterator(rhod_rv.begin(), ijk.begin()),
	    thrust::make_permutation_iterator(T.begin(),       ijk.begin()),
	    thrust::make_permutation_iterator(p.begin(),       ijk.begin()),
	    thrust::make_permutation_iterator(RH.begin(),      ijk.begin()),
	    rd3.begin(),
	    kpa.begin()
          )
        ), 
	rw2.begin(),                    // output
        detail::advance_rw2<real_t>(dt)
      );

      // calculating the 3rd wet moment after condensation (still not divided by dv)
      cond_dm3_helper();

      // adding the third moment after condensation to dm_3
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,            // input - 1st arg
        thrust::make_permutation_iterator(drhod_rv.begin(), count_ijk.begin()), // input - 2nd arg
        thrust::make_permutation_iterator(drhod_rv.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );

      // multiplying dv*dm_3 by -rho_w*4/3*pi/dv
      thrust::transform(
        drhod_rv.begin(), drhod_rv.end(),        // input - 1st arg
        thrust::make_constant_iterator<real_t>(  // input - 2nd arg
          - common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres
          * real_t(4./3) * pi<real_t>()
          / (opts.dx * opts.dy * opts.dz)
        ),
        drhod_rv.begin(),                        // output
        thrust::multiplies<real_t>()
      );

      // updating rhod_rv 
      thrust::transform(
        rhod_rv.begin(), rhod_rv.end(),  // input - 1st arg
        drhod_rv.begin(),                // input - 2nd arg
        rhod_rv.begin(),                 // output
        thrust::plus<real_t>() 
      );
      assert(*thrust::min_element(rhod_rv.begin(), rhod_rv.end()) >= 0);

      // updating rhod_th
      {
        typedef thrust::zip_iterator<thrust::tuple<
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator
        > > zip_it_t;
 
	thrust::transform(
	  rhod_th.begin(), rhod_th.end(),  // input - 1st arg
	  thrust::transform_iterator<      // input - 2nd arg
	    detail::drhod_th<real_t>,
	    zip_it_t,
	    real_t
	  >(
            zip_it_t(thrust::make_tuple(  // args
	      drhod_rv.begin(), // 
	      rhod.begin(),     // drhod_th = (drhod_rv / rhod) * d_rhodtheta_d_rv(T, rhod_th)
	      T.begin(),        //
	      rhod_th.begin()   //
	    )),
	    detail::drhod_th<real_t>()     // func
	  ),
	  rhod_th.begin(),                 // output
	  thrust::plus<real_t>()
	);
      }
    }
  };  
};
