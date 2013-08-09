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
        real_t operator()(const thrust::tuple<n_t, real_t> &tpl)
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
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t> &tpl)
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
      struct advance_rw2
      {
        quantity<si::time, real_t> dt;

        // ctor
        advance_rw2(const real_t &dt) : dt(dt * si::seconds) {}

	// rw2_new = rw2_old + f_rw2(rw2_new) * dt
	// rw2_new = rw2_old + 2 * rw * f_rw(rw2_new) * dt
// TODO: now it's Euler scheme!!
        BOOST_GPU_ENABLED
        real_t operator()(
          const real_t &rw2_unitless,
          const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t> &tpl
        ) {
          const quantity<si::length,        real_t> rw      = std::sqrt(rw2_unitless) * si::metres;
          const quantity<si::area,          real_t> rw2     = rw2_unitless * si::square_metres;
          const quantity<si::volume,        real_t> rw3     = rw * rw * rw;
          const quantity<si::mass_density,  real_t> rhod_rv = thrust::get<0>(tpl) * si::kilograms / si::cubic_metres;
          const quantity<si::temperature,   real_t> T       = thrust::get<1>(tpl) * si::kelvins;
	  const quantity<si::pressure,      real_t> p       = thrust::get<2>(tpl) * si::pascals;
	  const quantity<si::dimensionless, real_t> RH      = thrust::get<3>(tpl);
	  const quantity<si::volume,        real_t> rd3     = thrust::get<4>(tpl) * si::cubic_metres;
	  const quantity<si::dimensionless, real_t> kpa     = thrust::get<5>(tpl);

          using namespace common::maxwell_mason;
          using namespace common::kappa_koehler;
          using namespace common::kelvin;

          return (rw2 + dt * real_t(2) * rdrdt( 
            rhod_rv, T, p, RH, 
            a_w(rw3, rd3, kpa),
            klvntrm(rw, T)
          )) / si::square_metres;
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
	    rhod_rv.begin(),
	    T.begin(),
	    p.begin(),
	    RH.begin(),
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
