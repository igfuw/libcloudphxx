/** @file
  * @copyright University of Warsaw
  * @brief Autoconversion and collection righ-hand side terms using Kessler formulae
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
      cont_t &dot_ra_cont,
      const cont_t &rc_cont,
      const cont_t &rr_cont,
      const cont_t &ra_cont,
      const cont_t &th_cont,
      const cont_t &p_cont,
      const cont_t &rhod_cont,
      const real_t &dt
    )
//</listing>
    {
      for (auto tup : zip(dot_rc_cont, dot_rr_cont, dot_ra_cont, rc_cont, rr_cont, ra_cont, th_cont, p_cont, rhod_cont))
      {
        real_t
          rc_to_rr = 0,
          rc_to_ra = 0,
          &dot_rc = std::get<0>(tup),
          &dot_rr = std::get<1>(tup),
          &dot_ra = std::get<2>(tup);
        const real_t
          &rc     = std::get<3>(tup),
          &rr     = std::get<4>(tup),
          &ra     = std::get<5>(tup),
          &th     = std::get<6>(tup),
          &p      = std::get<7>(tup),
          &rhod   = std::get<8>(tup);

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

        // ice A heterogeneous nucleation
        //if (opts.hetA)
        //{
          rc_to_ra += (
            formulae::het_A_nucleation(
                    ra * si::dimensionless(),
                    rc * si::dimensionless(),
                    th * si::kelvin,
                	p * si::pascal,
                    rhod * si::mass_density,
                    dt * si::time
                  ) * si::seconds // to make it dimensionless
          );
        //}

        // ice A homogeneous nucleation 2
        //if (opts.homA2)
        //{
        rc_to_ra += (
          formulae::hom_A_nucleation_2(
                  rc * si::dimensionless(),
                  th * si::kelvin,
                  p * si::pascal,
                  dt * si::time
                ) * si::seconds // to make it dimensionless
        );
        //}

	dot_rr += rc_to_rr;
	dot_rc -= rc_to_rr;
  dot_rc -= rc_to_ra;
  dot_ra += rc_to_ra;
      }
    }

    template <typename real_t, class cont_t>
    void rhs_cellwise_nwtrph(
      const opts_t<real_t> &opts,
      cont_t &dot_th_cont,
      cont_t &dot_rv_cont,
      cont_t &dot_rc_cont,
      cont_t &dot_rr_cont,
      cont_t &dot_ra_cont,
      const cont_t &rhod_cont,
      const cont_t &p_cont,
      const cont_t &th_cont,
      const cont_t &rv_cont,
      const cont_t &rc_cont,
      const cont_t &rr_cont,
      const cont_t &ra_cont,
      const real_t &dt // time step in seconds
    )
    {
      rhs_cellwise<real_t, cont_t>(opts, dot_rc_cont, dot_rr_cont, dot_ra_cont, rc_cont, rr_cont, ra_cont, th_cont, p_cont, rhod_cont, dt); //calculate T from theta

      // rain evaporation treated as a force in Newthon-Raphson saturation adjustment
      for (auto tup : zip(dot_th_cont, dot_rv_cont, dot_rr_cont, rhod_cont, p_cont, th_cont, rv_cont, rr_cont))
      {
        using namespace common;

        real_t
          tmp = 0,
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

        tmp = (
          formulae::evaporation_rate(
            rv * si::dimensionless(),
            r_vs * si::dimensionless(),
            rr * si::dimensionless(),
            rhod,
            p
          ) * si::seconds * dt
        );

        // limiting
        tmp = std::min(rr, tmp) / dt;

        dot_rv += tmp;
        dot_rr -= tmp;

        dot_th -= const_cp::l_v(T) / (moist_air::c_pd<real_t>() * theta_std::exner(p)) * tmp / si::kelvins;
      }
    }
  }
}
