/** @file
  * @copyright University of Warsaw
  * @brief Autoconversion and collection righ-hand side terms using Kessler formulae.
  *        Ice nucleation, growth by deposition and riming, and melting from Grabowski (1999).
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include "extincl.hpp"

namespace libcloudphxx
{
  namespace blk_1m
  {
    //<listing>
    template <typename real_t, class cont_t>
    void rhs_cellwise(
      const opts_t<real_t> &opts,
      cont_t &dot_rc_cont,
      cont_t &dot_rr_cont,
      const cont_t &rc_cont,
      const cont_t &rr_cont
    )
//</listing>
    {
      for (auto tup : zip(dot_rc_cont, dot_rr_cont, rc_cont, rr_cont))
      {
        real_t
          rc_to_rr = 0,
          &dot_rc = std::get<0>(tup),
          &dot_rr = std::get<1>(tup);
        const real_t
          &rc     = std::get<2>(tup),
          &rr     = std::get<3>(tup);

        // autoconversion
        if (opts.conv)
        {
          rc_to_rr += (
            formulae::autoconversion_rate(
                    rc        * si::dimensionless(),
                    opts.r_c0 * si::dimensionless(),
                    opts.k_acnv / si::seconds
                  ) * si::seconds // to make it dimensionless
          );
        }

        // collection
        if (opts.accr)
        {
          rc_to_rr += (
            formulae::collection_rate(
                    rc * si::dimensionless(),
                    rr * si::dimensionless()
                  ) * si::seconds // to make it dimensionless
          );
        }

        dot_rr += rc_to_rr;
        dot_rc -= rc_to_rr;
      }
    }

    template <typename real_t, class cont_t>
    void rhs_cellwise_nwtrph(
      const opts_t<real_t> &opts,
      cont_t &dot_th_cont,
      cont_t &dot_rv_cont,
      cont_t &dot_rc_cont,
      cont_t &dot_rr_cont,
      const cont_t &rhod_cont,
      const cont_t &p_cont,
      const cont_t &th_cont,
      const cont_t &rv_cont,
      const cont_t &rc_cont,
      const cont_t &rr_cont,
      const real_t &dt
    )
    {
      rhs_cellwise<real_t, cont_t>(opts, dot_rc_cont, dot_rr_cont,  rc_cont, rr_cont);

      // rain evaporation treated as a force in Newthon-Raphson saturation adjustment
      for (auto tup : zip(dot_th_cont, dot_rv_cont, dot_rr_cont, rhod_cont, p_cont, th_cont, rv_cont, rr_cont))
      {
        using namespace common;

        real_t
          rr_to_rv = 0,
          &dot_th = std::get<0>(tup),
          &dot_rv = std::get<1>(tup),
          &dot_rr = std::get<2>(tup);

        const quantity<si::mass_density, real_t>
          rhod = std::get<3>(tup) * si::kilograms / si::cubic_metres;

        const quantity<si::pressure, real_t>
          p   = std::get<4>(tup) * si::pascals;

        const quantity<si::temperature, real_t>
          th = std::get<5>(tup) * si::kelvins;

        const real_t
          &rv = std::get<6>(tup),
          &rr = std::get<7>(tup);

        quantity<si::temperature, real_t> T = th * theta_std::exner(p);
        real_t r_vs = const_cp::r_vs(T, p);

        rr_to_rv += (
          formulae::evaporation_rate(
            rv * si::dimensionless(),
            r_vs * si::dimensionless(),
            rr * si::dimensionless(),
            rhod,
            p
          ) * si::seconds * dt
        );

        // limiting
        rr_to_rv = std::min(rr, rr_to_rv) / dt;

        dot_rv += rr_to_rv;
        dot_rr -= rr_to_rv;
        dot_th -= const_cp::l_v(T) / (moist_air::c_pd<real_t>() * theta_std::exner(p)) * dot_rv / si::kelvins;
      }
    }

    template <typename real_t, class cont_t>
    void rhs_cellwise_nwtrph_ice(
      const opts_t<real_t> &opts,
      cont_t &dot_th_cont,
      cont_t &dot_rv_cont,
      cont_t &dot_rc_cont,
      cont_t &dot_rr_cont,
      cont_t &dot_ria_cont,
      //cont_t &dot_rib_cont,
      const cont_t &rhod_cont,
      const cont_t &p_cont,
      const cont_t &th_cont,
      const cont_t &rv_cont,
      const cont_t &rc_cont,
      const cont_t &rr_cont,
      const cont_t &ria_cont,
      //const cont_t &rib_cont,
      const real_t &dt
    )
    {
      // autoconversion, collection and rain evaporation:
      rhs_cellwise_nwtrph<real_t, cont_t>(opts, dot_th_cont, dot_rv_cont,  dot_rc_cont, dot_rr_cont, rhod_cont, p_cont, th_cont, rv_cont, rc_cont, rr_cont, dt);

      // rain evaporation treated as a force in Newthon-Raphson saturation adjustment
      for (auto tup : zip(dot_th_cont, dot_rv_cont, dot_rc_cont, dot_rr_cont, dot_ria_cont, rhod_cont, p_cont, th_cont, rv_cont, rc_cont, rr_cont, ria_cont))
      {
        using namespace common;

        real_t
          rv_to_ria = 0,
          rc_to_ria = 0,
        &dot_th = std::get<0>(tup),
        &dot_rv = std::get<1>(tup),
        &dot_rc = std::get<2>(tup),
        &dot_rr = std::get<3>(tup),
        &dot_ria = std::get<4>(tup);

        const quantity<si::mass_density, real_t>
          rhod = std::get<5>(tup) * si::kilograms / si::cubic_metres;

        const quantity<si::pressure, real_t>
          p   = std::get<6>(tup) * si::pascals;

        const quantity<si::temperature, real_t>
          th = std::get<7>(tup) * si::kelvins;

        const real_t
          &rv     = std::get<8>(tup),
          &rc     = std::get<9>(tup),
          &rr     = std::get<10>(tup),
          &ria     = std::get<11>(tup);

        quantity<si::temperature, real_t> T = th * theta_std::exner(p);
        real_t rvs = const_cp::r_vs(T, p);
        real_t rvsi = const_cp::r_vsi(T, p);

        // ice A heterogeneous nucleation
        if (opts.hetA)
        {
          rc_to_ria += (
            formulae::het_A_nucleation(
                    ria * si::dimensionless(),
                    rc * si::dimensionless(),
                    T,
                    rhod,
                    dt * si::seconds
                  ) * si::seconds // to make it dimensionless
          );
        }

        // ice A homogeneous nucleation 1
        if (opts.homA1)
        {
          rv_to_ria += (
            formulae::hom_A_nucleation_1(
                    rv * si::dimensionless(),
                    rvs * si::dimensionless(),
                    rvsi * si::dimensionless(),
                    T,
                    dt * si::seconds
                  ) * si::seconds // to make it dimensionless
          );
        }

        // ice A homogeneous nucleation 2
        if (opts.homA2)
        {
          rc_to_ria += (
            formulae::hom_A_nucleation_2(
                    rc * si::dimensionless(),
                    T,
                    dt * si::seconds
                  ) * si::seconds // to make it dimensionless
          );
        }

        dot_rc -= rc_to_ria;
        dot_rv -= rv_to_ria;
        dot_ria += rc_to_ria + rv_to_ria;
        dot_th += const_cp::l_s(T) / (moist_air::c_pd<real_t>() * theta_std::exner(p)) * rv_to_ria / si::kelvins; //heat of sublimation
        dot_th += const_cp::l_f(T) / (moist_air::c_pd<real_t>() * theta_std::exner(p)) * rc_to_ria / si::kelvins; //heat of freezing
      }
    }

  }
}
