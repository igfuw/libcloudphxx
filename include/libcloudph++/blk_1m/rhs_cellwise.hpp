/** @file
  * @copyright University of Warsaw
  * @brief Autoconversion and collection righ-hand side terms using Kessler formulae
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/blk_1m/extincl.hpp>

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
          tmp = 0,
          &dot_rc = boost::get<0>(tup),
          &dot_rr = boost::get<1>(tup);
        const real_t
          &rc     = boost::get<2>(tup),
          &rr     = boost::get<3>(tup);

        // autoconversion
        if (opts.conv)
        {
	  tmp += ( 
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
	  tmp += (
	    formulae::collection_rate(
              rc * si::dimensionless(), 
              rr * si::dimensionless()
            ) * si::seconds // to make it dimensionless
	  );
        }

	dot_rr += tmp;
	dot_rc -= tmp;
      }
    }

    template <typename real_t, class cont_t>
    void rhs_cellwise_nwtrph(
      const opts_t<real_t> &opts,
      cont_t &dot_rc_cont, 
      cont_t &dot_rr_cont,
      const cont_t &rc_cont,
      const cont_t &rr_cont,
      const cont_t &rhod_cont,
      const cont_t &p_cont,
      const cont_t &p_d_cont,
      const cont_t &rv_cont,
      const cont_t &th_cont,
      cont_t &dot_rv_cont,
      cont_t &dot_th_cont
    )   
    {
      rhs_cellwise<real_t, cont_t>(opts, dot_rc_cont, dot_rr_cont, rc_cont, rr_cont);
      
      // rain evaporation treated as a force in Newthon-Raphson saturation adjustment
      for (auto tup : zip(dot_rr_cont, rr_cont, rhod_cont, p_cont, p_d_cont, rv_cont, th_cont, dot_rv_cont, dot_th_cont))
      {
        using namespace common;
        
        real_t
          tmp = 0,
          &dot_rr = boost::get<0>(tup);
        const real_t
          &rr     = boost::get<1>(tup);

        const quantity<si::mass_density, real_t> 
          rhod = boost::get<2>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::pressure, real_t> 
          p    = boost::get<3>(tup) * si::pascals,
          p_d  = boost::get<4>(tup) * si::pascals;
        const real_t
          rv   = boost::get<5>(tup);
        const quantity<si::temperature, real_t> 
          th   = boost::get<6>(tup) * si::kelvins;
        real_t
          &dot_rv = boost::get<7>(tup),
          &dot_th = boost::get<8>(tup);

        quantity<si::temperature, real_t> T = th * theta_std::exner(p_d);
        real_t r_vs = const_cp::r_vs(T, p_d);

        tmp = (
          formulae::evaporation_rate(
            rv * si::dimensionless(), 
            r_vs * si::dimensionless(),
            rr * si::dimensionless(),
            rhod,
            p
          ) * si::seconds // to make it dimensionless
        );

        dot_rv += tmp;
        dot_rr -= tmp;

        dot_th -= const_cp::l_v(T) / (moist_air::c_pd<real_t>() * theta_std::exner(p_d)) * tmp / si::kelvins;
      }
    }
  }
}
