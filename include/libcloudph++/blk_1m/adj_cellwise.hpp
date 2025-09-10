/** @file
  * @copyright University of Warsaw
  * @brief saturation adjustment routine using Boost.odeint for
  *        solving latent-heat release equation
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include "extincl.hpp"

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
        bool const_p, // true if a constant pressure profile is used (e.g. anelastic model)
             th_dry; // true if th is dry potential temperature; "standard" potential temperature otherwise

        public:

        quantity<si::pressure,     real_t> p;   // total pressure

        void init(
          const quantity<si::mass_density, real_t> &_rhod,
          const quantity<si::pressure, real_t> &_p,
          const quantity<si::temperature, real_t> &th,
          const quantity<si::dimensionless, real_t> &rv,
          const bool _const_p,
          const bool _th_dry
        )
        {
          const_p = _const_p;
          th_dry = _th_dry;
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

          if(!const_p && th_dry)
          {
            T = common::theta_dry::T<real_t>(th, rhod);
            p = common::theta_dry::p<real_t>(rhod, rv, T);
          }
          else if(const_p && !th_dry)
          {
            T = th * common::theta_std::exner(p);
          }
          else throw std::runtime_error("adj_cellwise: one (and only one) of opts.const_p and opts.th_dry must be true");

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
      const cont_t &rhod_cont,
      const cont_t &p_cont,
      cont_t &th_cont,
      cont_t &rv_cont,
      cont_t &rc_cont,
      const real_t &dt
    )
    {
      using namespace common;
      if (!opts.cond) return; // ignoring values of opts.cevp

      for (auto tup : zip(rhod_cont, p_cont, th_cont, rv_cont, rc_cont))
      {
        real_t
          &th = std::get<2>(tup),
          &rv = std::get<3>(tup),
          &rc = std::get<4>(tup);

        // double-checking....
        //assert(th >= 273.15); // TODO: that's theta, not T!
        assert(rc >= 0);
        assert(rv >= 0);

        real_t drc = 0;

        real_t rv_tmp = rv;
        quantity<si::temperature, real_t> th_tmp = th * si::kelvins;

        quantity<si::temperature, real_t>   T, T_tmp;
        quantity<si::mass_density, real_t>  rhod;
        quantity<si::pressure, real_t>      p;
        quantity<si::dimensionless, real_t> exner;

        if(!opts.const_p && opts.th_dry)
        {
          rhod = std::get<2>(tup) * si::kilograms / si::cubic_metres;
          T = common::theta_dry::T<real_t>(th_tmp, rhod);
          p = common::theta_dry::p<real_t>(rhod, rv, T);
        }
        else if(opts.const_p && !opts.th_dry)
        {
          p = std::get<1>(tup) * si::pascals;
          exner = common::theta_std::exner(p);
          T = th_tmp * exner;
        }
        else throw std::runtime_error("adj_cellwise: one (and only one) of opts.const_p and opts.th_dry must be true");

        // constant l_v used in theta update
        auto L0 = const_cp::l_v(T);

        T_tmp = T;
        for (int iter = 0; iter < opts.nwtrph_iters; ++iter)
        {
          // TODO: use the approximate Tetens formulas for p_vs and r_vs from tetens.hpp?
          quantity<si::pressure, real_t> p_vs = const_cp::p_vs(T_tmp);

          // tricky, constant L0 comes from theta = theta + L0 / (c_pd * exner_p) * drc
          // while variable L comes from dp_vs/dT
          auto L = const_cp::l_v(T_tmp);
          real_t coeff = L * L0 / (moist_air::c_pd<real_t>() * moist_air::R_v<real_t>()) / (T_tmp * T_tmp) / (1 - p_vs / p);

          real_t r_vs = const_cp::r_vs(T_tmp, p);

          drc +=  (rv_tmp - r_vs) / (1 + coeff * r_vs);

          rv_tmp = rv - drc;
          th_tmp = th * si::kelvins + th_tmp / T_tmp * L0 / (moist_air::c_pd<real_t>()) * drc;

          if(opts.th_dry)
            T_tmp = common::theta_dry::T<real_t>(th_tmp, rhod);
          else
            T_tmp = th_tmp * exner;
        }

        // limiting
        drc = std::min(rv, std::max(-rc, drc));

        rv -= drc;
        rc += drc;
        th += th / T * L0 / (moist_air::c_pd<real_t>()) * drc;

        // triple-checking....
        //assert(th >= 273.15); // that is theta, not T ! TODO
        assert(rc >= 0);
        assert(rv >= 0);
      }
    }
/*
// saturation adjustment with variable pressure calculated using rhod, rv and T; and th=th_dry
//<listing>
    template <typename real_t, class cont_t>
    void adj_cellwise_nwtrph(
      const opts_t<real_t> &opts,
      const cont_t &p_cont,
      cont_t &th_cont,
      cont_t &rv_cont,
      cont_t &rc_cont,
      const real_t &dt
    )
//</listing>
    {
      adj_cellwise_nwtrph_hlpr(opts, p_cont, th_cont, rv_cont, rc_cont, rr_cont, dt, false);
    }

// saturation adjustment with a constant pressure profile (e.g. anleastic model) and th=th_std
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
*/

// Below are saturation adjustment functions that use RK4 integration (from boost.odeint)

//<listing>
    template <typename real_t, class cont_t>
    void adj_cellwise_rk4(
      const opts_t<real_t> &opts,
      const cont_t &rhod_cont,
      const cont_t &p_cont, // value not used if opts.const_p = false
      cont_t &th_cont, // if opts._thdry, its th_dry, otherwise its th_std
      cont_t &rv_cont,
      cont_t &rc_cont,
      cont_t &rr_cont,
      const real_t &dt
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
          &rhod = std::get<0>(tup),
          &p    = std::get<1>(tup);
        real_t
          &th = std::get<2>(tup),
          &rv = std::get<3>(tup),
          &rc = std::get<4>(tup),
          &rr = std::get<5>(tup);

        // double-checking....
        //assert(th >= 273.15); // TODO: that's theta, not T!
        assert(rc >= 0);
        assert(rv >= 0);
        assert(rr >= 0);

        F.init(
          rhod * si::kilograms / si::cubic_metres,
          p    * si::pascals,
          th   * si::kelvins,
          rv   * si::dimensionless(),
          opts.const_p,
          opts.th_dry
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
        //assert(th >= 273.15); // that is theta, not T ! TODO
        assert(rc >= 0);
        assert(rv >= 0);
        assert(rr >= 0);
      }
    }


//<listing>
    template <typename real_t, class cont_t>
    void adj_cellwise(
      const opts_t<real_t> &opts,
      const cont_t &rhod_cont, // used only if opts.th_dry == true
      const cont_t &p_cont, // used only if opts.const_p == true
      cont_t &th_cont,
      cont_t &rv_cont,
      cont_t &rc_cont,
      cont_t &rr_cont, // used only if opts_adj_nwtrph == false
      const real_t &dt
    )
//</listing>
    {
      if(opts.adj_nwtrph)
        adj_cellwise_nwtrph(opts, rhod_cont, p_cont, th_cont, rv_cont, rc_cont, dt);
      else
        adj_cellwise_rk4(opts, rhod_cont, p_cont, th_cont, rv_cont, rc_cont, rr_cont, dt);
    }
  }

  /*

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
      */
};
