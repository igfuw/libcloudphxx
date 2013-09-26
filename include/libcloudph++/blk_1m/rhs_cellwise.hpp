/** @file
  * @copyright University of Warsaw
  * @brief Autoconversion and collection righ-hand side terms using Kessler formulae
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <algorithm>
#include <libcloudph++/blk_1m/formulae.hpp>

namespace libcloudphxx
{
  namespace blk_1m
  {
//<listing>
    template <typename real_t, class cont_t>
    void rhs_cellwise(
      const opts_t<real_t> &opts,
      cont_t &dot_rho_c_cont, 
      cont_t &dot_rho_r_cont,
      const cont_t &rho_d_cont,   
      const cont_t &rho_c_cont,
      const cont_t &rho_r_cont
    )   
//</listing>
    {
      for (auto tup : zip(dot_rho_c_cont, dot_rho_r_cont, rho_d_cont, rho_c_cont, rho_r_cont))
      {
        real_t
          tmp = 0,
          &dot_rho_c = boost::get<0>(tup),
          &dot_rho_r = boost::get<1>(tup);
        const real_t
          &rho_d     = boost::get<2>(tup),
          &rho_c  = boost::get<3>(tup),
          &rho_r  = boost::get<4>(tup);

        // autoconversion
        if (opts.conv)
        {
	  tmp += rho_d * ( 
	    formulae::autoconversion_rate(rho_c / rho_d * si::dimensionless(), opts.r_c0 * si::dimensionless()) 
	    * si::seconds // to make it dimensionless
	  );
        }

        // collection
        if (opts.accr)
        {
	  tmp += rho_d * (
	    formulae::collection_rate(rho_c / rho_d * si::dimensionless(), rho_r/rho_d * si::dimensionless()) * si::seconds
	  );
        }

	dot_rho_r += tmp;
	dot_rho_c -= tmp;
      }
    }    
  };
};
