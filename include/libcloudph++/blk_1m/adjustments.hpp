/** @file
  * @copyright University of Warsaw
  * @brief saturation adjustment routine using Boost.odeint for
  *        solving latent-heat release equation
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <boost/numeric/odeint.hpp>

#include <libcloudph++/common/const_cp.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/detail/zip.hpp>

#include <libcloudph++/blk_1m/formulae.hpp>

namespace libcloudphxx
{
  namespace blk_1m
  {
    namespace detail
    {
      // ODE rhs describing latent-heat release
      template <typename real_t>
      class rhs
      {
	private: 
         
        quantity<si::mass_density, real_t> rhod;
	
        public: 
        
        void init(
	  const quantity<si::mass_density, real_t> _rhod,
	  const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th,
	  const quantity<si::mass_density, real_t> rhod_rv
	)
	{
	  rhod = _rhod;
	  update(rhod_th, rhod_rv);
	}

	quantity<si::dimensionless, real_t> r, rs;
	quantity<si::pressure,      real_t> p;
	quantity<si::temperature,   real_t> T;

	private: 

        void update(
	  const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th,
	  const quantity<si::mass_density, real_t> rhod_rv
	)
	{
	  r  = rhod_rv / rhod;
	  T  = common::theta_dry::T<real_t>(rhod_th, rhod);
	  p  = common::theta_dry::p<real_t>(rhod, r, T);
	  rs = common::const_cp::r_vs<real_t>(T, p);
	}

	public: 

	// F = d (rho_d * th) / d (rho_d * r) 
        void operator()(
	  const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th,
	  quantity<si::temperature, real_t> &F,
	  const quantity<si::mass_density, real_t> rhod_rv
	)
	{
	  update(rhod_th, rhod_rv);
	  F = common::theta_dry::d_rhodtheta_d_rv<real_t>(T, rhod_th) / rhod; 
	}
      };
    }    

//<listing>
    template <typename real_t, class cont_t>
    void adjustments(
      const opts_t<real_t> &opt,
      const cont_t &rhod_cont, 
      cont_t &rhod_th_cont, 
      cont_t &rhod_rv_cont,
      cont_t &rhod_rc_cont,
      cont_t &rhod_rr_cont
    )
//</listing>
    {
      assert(opt.dt != 0);

      namespace odeint = boost::numeric::odeint;

      // odeint::euler< // TODO: opcja?
      odeint::runge_kutta4<
	quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>, // state_type
	real_t, // value_type
	quantity<si::temperature, real_t>, // deriv_type
	quantity<si::mass_density, real_t>, // time_type
	odeint::vector_space_algebra,
	odeint::default_operations,
	odeint::never_resizer
      > S; // TODO: would be better to instantiate in the ctor (but what about thread safety! :()
      typename detail::rhs<real_t> F;

      for (auto tup : zip(rhod_cont, rhod_th_cont, rhod_rv_cont, rhod_rc_cont, rhod_rr_cont))
      {
        const real_t
          &rhod = boost::get<0>(tup);
        real_t 
          &rhod_th = boost::get<1>(tup), 
          &rhod_rv = boost::get<2>(tup), 
          &rhod_rc = boost::get<3>(tup), 
          &rhod_rr = boost::get<4>(tup);

	F.init(
	  rhod    * si::kilograms / si::cubic_metres,
	  rhod_th * si::kilograms / si::cubic_metres * si::kelvins,
	  rhod_rv * si::kilograms / si::cubic_metres
	);

	real_t // TODO: quantity<si::mass_density
	  rho_eps = .00002, // TODO: as an option? // TODO: is it OK that we use rho_eps and integrate over r_v?
	  vapour_excess;
	real_t drho_rr_max = 0; // TODO: quantity<si::mass_density
	if (F.rs > F.r && rhod_rr > 0 && opt.revp)
	  drho_rr_max = 
            (opt.dt * si::seconds) * formulae::evaporation_rate(F.r, F.rs, rhod_rr * si::kilograms / si::cubic_metres, F.p)
            / si::kilograms * si::cubic_metres; // to make it dimensionless
	bool incloud;

	// TODO: rethink and document rho_eps!!!
        // TODO: rho_eps as parameter
	while (
	  // condensation of cloud water if supersaturated more than threshold
	  (vapour_excess = rhod_rv - rhod * F.rs) > rho_eps
	  || 
          ( 
            opt.cevp && vapour_excess < -rho_eps && ( // or if subsaturated and 
	      (incloud = (rhod_rc > 0)) // in cloud (then cloud evaporation first)
	      ||                        // or 
              (opt.revp && rhod_rr > 0) // in rain shaft (rain evaporation out-of-cloud)
	    )
          )
	)
	{
          // initial guess for drho_rv
	  real_t drho_rv = - copysign(.5 * rho_eps, vapour_excess); // TODO: .5 - arbitrary!!! 
          // preventing negative mixing ratios if evaporating
	  if (vapour_excess < 0) drho_rv = 
            incloud ? std::min(rhod_rc, drho_rv) // limiting by rhod_rc
	            : std::min(drho_rr_max, std::min(rhod_rr, drho_rv)); // limiting by rhod_rr and drho_rr_max
	  assert(drho_rv != 0); // otherwise it should not pass the while condition!

	  // theta is modified by do_step, and hence we cannot pass an expression and we need a temp. var.
	  quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>
	    tmp = rhod_th * si::kilograms / si::cubic_metres * si::kelvins;

	  // integrating the First Law for moist air
	  S.do_step(
	    boost::ref(F),
	    tmp,
	    rhod_rv * si::kilograms / si::cubic_metres,
	    drho_rv * si::kilograms / si::cubic_metres
	  );

	  // latent heat source/sink due to evaporation/condensation
	  rhod_th = tmp / (si::kilograms / si::cubic_metres * si::kelvins);

	  // updating rhod_rv
	  rhod_rv += drho_rv;
	  assert(rhod_rv >= 0);
	  
	  if (vapour_excess > 0 || incloud)
	  {
            // condensation or evaporation of cloud water
	    rhod_rc -= drho_rv;
	    assert(rhod_rc >= 0);
	  }
	  else 
	  {
            // evaporation of rain water
	    assert(opt.revp); // should be guaranteed by the while() condition above
	    rhod_rr -= drho_rv;
	    assert(rhod_rr >= 0);
	    if ((drho_rr_max -= drho_rv) == 0) break; // but not more than Kessler allows
	  }
	}

	// hopefully true for RK4
	assert(F.r == real_t(rhod_rv / rhod));
	// double-checking....
	assert(rhod_rc >= 0);
	assert(rhod_rv >= 0);
	assert(rhod_rr >= 0);
      }
    }
  }
};
