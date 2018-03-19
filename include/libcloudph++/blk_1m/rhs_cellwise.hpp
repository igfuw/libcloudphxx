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
  };
};
