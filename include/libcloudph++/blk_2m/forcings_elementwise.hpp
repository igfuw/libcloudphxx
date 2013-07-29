/** @file
  * @copyright University of Warsaw
  * @brief Autoconversion and collection righ-hand side terms using Kessler formulae
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <algorithm>
#include <libcloudph++/blk_2m/formulae.hpp>

namespace libcloudphxx
{
  namespace blk_2m
  {
    template <typename real_t, class container_t>
    void forcings_elementwise(
      const opts_t<real_t> &opt,
      container_t drhod_th_cont,
      container_t drhod_rv_cont,
      container_t drhod_rc_cont,
      container_t drhod_rr_cont,
      container_t drhod_nc_cont,
      container_t drhod_nr_cont,
      const container_t rhod_cont,   
      const container_t rhod_th_cont,
      const container_t rhod_rv_cont,
      const container_t rhod_rc_cont,
      const container_t rhod_rr_cont,
      const container_t rhod_nc_cont,
      const container_t rhod_nr_cont
    )   
    {
      for (auto tup : zip(drhod_rc_cont, drhod_rr_cont, rhod_cont, rhod_rc_cont, rhod_rr_cont))
      {
        real_t
          tmp = 0,
          &drhod_rc = boost::get<0>(tup),
          &drhod_rr = boost::get<1>(tup);
        const real_t
          &rhod     = boost::get<2>(tup),
          &rhod_rc  = boost::get<3>(tup),
          &rhod_rr  = boost::get<4>(tup);

        // autoconversion
        if (opt.conv)
        { 
	  tmp += rhod * ( 
	    formulae::autoconversion_rate(rhod_rc / rhod * si::dimensionless()) 
	    * si::seconds // to make it dimensionless
	  );
        }

        // collection
        if (opt.clct)
        {
	  tmp += rhod * (
	    formulae::collection_rate(rhod_rc/rhod * si::dimensionless(), rhod_rr/rhod * si::dimensionless()) * si::seconds
	  );
        }

	drhod_rr += tmp;
	drhod_rc -= tmp;
      }
    }    
  }
};
