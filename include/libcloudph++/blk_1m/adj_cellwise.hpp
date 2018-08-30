/** @file
  * @copyright University of Warsaw
  * @brief saturation adjustment routine using Boost.odeint for
  *        solving latent-heat release equation
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/blk_1m/extincl.hpp>

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
        bool const_p; // true if a constant pressure profile is used (e.g. anelastic model), otherwise pressure is deduced from rhod, rv and T
	
        public: 

	quantity<si::pressure,     real_t> p;   // total pressure
        
        void init(
	  const quantity<si::mass_density, real_t> &_rhod,
	  const quantity<si::pressure, real_t> &_p,
	  const quantity<si::temperature, real_t> &th,
	  const quantity<si::dimensionless, real_t> &rv,
          const bool _const_p
	)
	{
          const_p = _const_p;
	  rhod = _rhod;
          p    = _p;
	  update(th, rv);
	}

	quantity<si::dimensionless, real_t> r, rs;
	quantity<si::temperature,   real_t> T;

	private: 

        void update(
	  const quantity<si::temperature, real_t> &th,
	  const quantity<si::dimensionless, real_t> &rv
	)
	{
          r  = rv;
          if(!const_p)
          {
	    T  = common::theta_dry::T<real_t>(th, rhod);
	    p   = common::theta_dry::p<real_t>(rhod, rv, T);
          }
          else
          {
	    T  = th * common::theta_std::exner<real_t>(p);
          }

	  rs = common::const_cp::r_vs<real_t>(T, p);
	}

	public: 

	// F = d th / d rv 
        void operator()(
	  const quantity<si::temperature, real_t> &th,
	  quantity<si::temperature, real_t> &F,
	  const quantity<si::dimensionless, real_t> &rv
	)
	{
	  update(th, rv);
	  F = common::theta_dry::d_th_d_rv<real_t>(T, th); 
	}
      };
    }    

    template <typename real_t, class cont_t>
    void adj_cellwise_nwtrph(
      const opts_t<real_t> &opts,
      const cont_t &p_cont,
      cont_t &th_cont, 
      cont_t &rv_cont,
      cont_t &rc_cont,
      const real_t &dt
    )
    {
      using namespace common;
      if (!opts.cond) return; // ignoring values of opts.cevp

      for (auto tup : zip(p_cont, th_cont, rv_cont, rc_cont))
      {
        const quantity<si::pressure, real_t> 
          p    = boost::get<0>(tup) * si::pascals;

        real_t 
          &th = boost::get<1>(tup), 
          &rv = boost::get<2>(tup), 
          &rc = boost::get<3>(tup);
	
        // double-checking....
	assert(th >= 273.15); // TODO: that's theta, not T!
	assert(rc >= 0);
	assert(rv >= 0);

        real_t drc = 0;

        real_t rv_tmp = rv;
        quantity<si::temperature, real_t> th_tmp = th * si::kelvins;
	
        auto exner_p = theta_std::exner(p);
        
        quantity<si::temperature, real_t> T = th_tmp * exner_p; 

        // constant l_v used in theta update
        auto L0 = const_cp::l_v(T);

        for (int iter = 0; iter < opts.nwtrph_iters; ++iter)
        {
          // TODO: use the approximate Tetens formulas for p_vs and r_vs from tetens.hpp?
	  quantity<si::pressure, real_t> p_vs = const_cp::p_vs(T); 

          // tricky, constant L0 comes from theta = theta + L0 / (c_pd * exner_p) * drc
          // while variable L comes from dp_vs/dT
          auto L = const_cp::l_v(T);
          real_t coeff = L * L0 / (moist_air::c_pd<real_t>() * moist_air::R_v<real_t>()) / (T * T) / (1 - p_vs / p);

          real_t r_vs = const_cp::r_vs(T, p);

          drc +=  (rv_tmp - r_vs) / (1 + coeff * r_vs);

          rv_tmp = rv - drc;
          th_tmp = th * si::kelvins + L0 / (moist_air::c_pd<real_t>() * exner_p) * drc;
	  
          T = th_tmp * exner_p; 
        }

        // limiting
        drc = std::min(rv, std::max(-rc, drc));

        rv -= drc;
        rc += drc;
        th += L0 / (moist_air::c_pd<real_t>() * exner_p) * drc / si::kelvins;

	// triple-checking....
	assert(th >= 273.15); // that is theta, not T ! TODO
	assert(rc >= 0);
	assert(rv >= 0);
      }
    }

//<listing>
    template <typename real_t, class cont_t>
    void adj_cellwise_hlpr(
      const opts_t<real_t> &opts,
      const cont_t &rhod_cont, 
      const cont_t &p_cont, // value not used if const_p = false
      cont_t &th_cont, // if const_p == false, its th_dry, otherwise its th_std 
      cont_t &rv_cont,
      cont_t &rc_cont,
      cont_t &rr_cont,
      const real_t &dt,
      const bool const_p
    )
//</listing>
    {
      if (!opts.cond) return; // ignoring values of opts.cevp and opts.revp

      namespace odeint = boost::numeric::odeint;

      // odeint::euler< // TODO: opcja?
      odeint::runge_kutta4<
	quantity<si::temperature, real_t>,   // state_type
	real_t,                              // value_type
	quantity<si::temperature, real_t>,   // deriv_type
	quantity<si::dimensionless, real_t>, // time_type
	odeint::vector_space_algebra,
	odeint::default_operations,
	odeint::never_resizer
      > S; // TODO: would be better to instantiate in the ctor (but what about thread safety! :()
      typename detail::rhs<real_t> F;

      for (auto tup : zip(rhod_cont, p_cont, th_cont, rv_cont, rc_cont, rr_cont))
      {
        const real_t
          &rhod = boost::get<0>(tup),
          &p    = boost::get<1>(tup);
        real_t 
          &th = boost::get<2>(tup), 
          &rv = boost::get<3>(tup), 
          &rc = boost::get<4>(tup), 
          &rr = boost::get<5>(tup);

	// double-checking....
	assert(th >= 273.15); // TODO: that's theta, not T!
	assert(rc >= 0);
	assert(rv >= 0);
	assert(rr >= 0); 

	F.init(
	  rhod * si::kilograms / si::cubic_metres,
          p    * si::pascals,
	  th   * si::kelvins,
	  rv   * si::dimensionless(),
          const_p
	);

	real_t vapour_excess;
	real_t drr_max = 0;
	if (F.rs > F.r && rr > 0 && opts.revp) 
        {
          drr_max = (dt * si::seconds) * formulae::evaporation_rate(
            F.r, F.rs, rr * si::dimensionless(), rhod * si::kilograms / si::cubic_metres, F.p
	  );
        }
	bool incloud;

	// TODO: rethink and document r_eps!!!
	while (
	  // condensation of cloud water if supersaturated more than a threshold
	  (vapour_excess = rv - F.rs) > opts.r_eps
	  || 
          ( 
            opts.cevp && vapour_excess < -opts.r_eps && ( // or if subsaturated and 
	      (incloud = (rc > 0))  // in cloud (then cloud evaporation first)
	      ||                    // or 
              (opts.revp && rr > 0 && drr_max > 0) // in rain shaft (rain evaporation out-of-cloud)
	    )
          )
	)
	{
          // an arbitrary initial guess for drv
	  real_t drv = - copysign(std::min(.5 * opts.r_eps, .5 * vapour_excess), vapour_excess); 
          // preventing negative mixing ratios if evaporating
	  if (vapour_excess < 0) drv = 
            incloud ? std::min(rc, drv) // limiting by rc
	            : std::min(drr_max, std::min(rr, drv)); // limiting by rr and drr_max
	  assert(drv != 0); // otherwise it should not pass the while condition!

	  // theta is modified by do_step, and hence we cannot pass an expression and we need a temp. var.
	  quantity<si::temperature, real_t> tmp = th * si::kelvins;

	  // integrating the First Law for moist air
	  S.do_step(
	    boost::ref(F),
	    tmp,
	    rv  * si::dimensionless(),
	    drv * si::dimensionless()
	  );

	  // latent heat source/sink due to evaporation/condensation
	  th = tmp / si::kelvins;

	  // updating rv
	  rv += drv;
	  assert(rv >= 0);
	  
	  if (vapour_excess > 0 || incloud)
	  {
            // condensation or evaporation of cloud water
	    rc -= drv;
	    assert(rc >= 0);
	  }
	  else 
	  {
            // evaporation of rain water
	    assert(opts.revp); // should be guaranteed by the while() condition above
	    rr -= drv;
	    assert(rr >= 0);
	    if ((drr_max -= drv) == 0) break; // but not more than Kessler allows
	  }
	}
      
	// hopefully true for RK4
	assert(F.r == rv);
	// triple-checking....
	assert(th >= 273.15); // that is theta, not T ! TODO
	assert(rc >= 0);
	assert(rv >= 0);
	assert(rr >= 0);
      }
    }

// saturation adjustment with variable pressure calculated using rhod, rv and T
//<listing>
    template <typename real_t, class cont_t>
    void adj_cellwise(
      const opts_t<real_t> &opts,
      const cont_t &rhod_cont, 
      cont_t &th_cont,  // th_dry
      cont_t &rv_cont,
      cont_t &rc_cont,
      cont_t &rr_cont,
      const real_t &dt
    )
//</listing>
    {
      // rhod_cont passed twice on purpose - it's a dummy var for pressures to make the tuple compile
      adj_cellwise_hlpr(opts, rhod_cont, rhod_cont, th_cont, rv_cont, rc_cont, rr_cont, dt, false);
    }

// saturation adjustment with a constant pressure profile (e.g. anleastic model)
// needs a different name, because boost python got confused (TODO: fix it)
//<listing>
    template <typename real_t, class cont_t>
    void adj_cellwise_constp(
      const opts_t<real_t> &opts,
      const cont_t &rhod_cont, 
      const cont_t &p_cont, 
      cont_t &th_cont,  // th_std
      cont_t &rv_cont,
      cont_t &rc_cont,
      cont_t &rr_cont,
      const real_t &dt
    )
//</listing>
    {
      adj_cellwise_hlpr(opts, rhod_cont, p_cont, th_cont, rv_cont, rc_cont, rr_cont, dt, true);
    }
  }
};
