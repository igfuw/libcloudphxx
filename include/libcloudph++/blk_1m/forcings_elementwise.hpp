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
    void forcings_elementwise(
      const opts_t<real_t> &opt,
      cont_t &dot_rhod_rc_cont, 
      cont_t &dot_rhod_rr_cont,
      const cont_t &rhod_cont,   
      const cont_t &rhod_rc_cont,
      const cont_t &rhod_rr_cont
    )   
//</listing>
    {
      for (auto tup : zip(dot_rhod_rc_cont, dot_rhod_rr_cont, rhod_cont, rhod_rc_cont, rhod_rr_cont))
      {
        real_t
          tmp = 0,
          &dot_rhod_rc = boost::get<0>(tup),
          &dot_rhod_rr = boost::get<1>(tup);
        const real_t
          &rhod     = boost::get<2>(tup),
          &rhod_rc  = boost::get<3>(tup),
          &rhod_rr  = boost::get<4>(tup);

        // autoconversion
        if (opt.conv)
        {
	  tmp += rhod * ( 
	    formulae::autoconversion_rate(rhod_rc / rhod * si::dimensionless(), opt.r_c0 * si::dimensionless()) 
	    * si::seconds // to make it dimensionless
	  );
        }

        // collection
        if (opt.accr)
        {
	  tmp += rhod * (
	    formulae::collection_rate(rhod_rc/rhod * si::dimensionless(), rhod_rr/rhod * si::dimensionless()) * si::seconds
	  );
        }

	dot_rhod_rr += tmp;
	dot_rhod_rc -= tmp;
      }
    }    
  };
};
