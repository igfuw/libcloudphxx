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
         
        quantity<si::mass_density, real_t> rho_d;
	
        public: 
        
        void init(
	  const quantity<si::mass_density, real_t> _rho_d,
	  const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rho_e,
	  const quantity<si::mass_density, real_t> rho_v
	)
	{
	  rho_d = _rho_d;
	  update(rho_e, rho_v);
	}

	quantity<si::dimensionless, real_t> r, rs;
	quantity<si::pressure,      real_t> p;
	quantity<si::temperature,   real_t> T;

	private: 

        void update(
	  const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rho_e,
	  const quantity<si::mass_density, real_t> rho_v
	)
	{
	  r  = rho_v / rho_d;
	  T  = common::theta_dry::T<real_t>(rho_e, rho_d);
	  p  = common::theta_dry::p<real_t>(rho_d, r, T);
	  rs = common::const_cp::r_vs<real_t>(T, p);
	}

	public: 

	// F = d (rho_d * th) / d (rho_d * r) 
        void operator()(
	  const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rho_e,
	  quantity<si::temperature, real_t> &F,
	  const quantity<si::mass_density, real_t> rho_v
	)
	{
	  update(rho_e, rho_v);
	  F = common::theta_dry::d_rhodtheta_d_rv<real_t>(T, rho_e) / rho_d; 
	}
      };
    }    

//<listing>
    template <typename real_t, class cont_t>
    void adj_cellwise(
      const opts_t<real_t> &opts,
      const cont_t &rho_d_cont, 
      cont_t &rho_e_cont, 
      cont_t &rho_v_cont,
      cont_t &rho_c_cont,
      cont_t &rho_r_cont,
      const real_t &dt
    )
//</listing>
    {
      if (!opts.cond) return; // ignoring values of opts.cevp and opts.revp

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

      for (auto tup : zip(rho_d_cont, rho_e_cont, rho_v_cont, rho_c_cont, rho_r_cont))
      {
        const real_t
          &rho_d = boost::get<0>(tup);
        real_t 
          &rho_e = boost::get<1>(tup), 
          &rho_v = boost::get<2>(tup), 
          &rho_c = boost::get<3>(tup), 
          &rho_r = boost::get<4>(tup);

	F.init(
	  rho_d  * si::kilograms / si::cubic_metres,
	  rho_e  * si::kilograms / si::cubic_metres * si::kelvins,
	  rho_v  * si::kilograms / si::cubic_metres
	);

	real_t vapour_excess;
	real_t drho_rr_max = 0; // TODO: quantity<si::mass_density
	if (F.rs > F.r && rho_r > 0 && opts.revp)
	  drho_rr_max = 
            (dt * si::seconds) * formulae::evaporation_rate(F.r, F.rs, rho_r * si::kilograms / si::cubic_metres, F.p)
            / si::kilograms * si::cubic_metres; // to make it dimensionless
	bool incloud;

	// TODO: rethink and document rho_eps!!!
	while (
	  // condensation of cloud water if supersaturated more than threshold
	  (vapour_excess = rho_v - rho_d * F.rs) > opts.rho_eps
	  || 
          ( 
            opts.cevp && vapour_excess < -opts.rho_eps && ( // or if subsaturated and 
	      (incloud = (rho_c > 0)) // in cloud (then cloud evaporation first)
	      ||                        // or 
              (opts.revp && rho_r > 0) // in rain shaft (rain evaporation out-of-cloud)
	    )
          )
	)
	{
          // initial guess for drho_v
	  real_t drho_v = - copysign(.5 * opts.rho_eps, vapour_excess); // TODO: .5 - arbitrary!!! 
          // preventing negative mixing ratios if evaporating
	  if (vapour_excess < 0) drho_v = 
            incloud ? std::min(rho_c, drho_v) // limiting by rho_c
	            : std::min(drho_rr_max, std::min(rho_r, drho_v)); // limiting by rho_r and drho_rr_max
	  assert(drho_v != 0); // otherwise it should not pass the while condition!

	  // theta is modified by do_step, and hence we cannot pass an expression and we need a temp. var.
	  quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>
	    tmp = rho_e * si::kilograms / si::cubic_metres * si::kelvins;

	  // integrating the First Law for moist air
	  S.do_step(
	    boost::ref(F),
	    tmp,
	    rho_v * si::kilograms / si::cubic_metres,
	    drho_v * si::kilograms / si::cubic_metres
	  );

	  // latent heat source/sink due to evaporation/condensation
	  rho_e = tmp / (si::kilograms / si::cubic_metres * si::kelvins);

	  // updating rho_v
	  rho_v += drho_v;
	  assert(rho_v >= 0);
	  
	  if (vapour_excess > 0 || incloud)
	  {
            // condensation or evaporation of cloud water
	    rho_c -= drho_v;
	    assert(rho_c >= 0);
	  }
	  else 
	  {
            // evaporation of rain water
	    assert(opts.revp); // should be guaranteed by the while() condition above
	    rho_r -= drho_v;
	    assert(rho_r >= 0);
	    if ((drho_rr_max -= drho_v) == 0) break; // but not more than Kessler allows
	  }
	}

	// hopefully true for RK4
	assert(F.r == real_t(rho_v / rho_d));
	// double-checking....
	assert(rho_c >= 0);
	assert(rho_v >= 0);
	assert(rho_r >= 0);
      }
    }
  }
};
