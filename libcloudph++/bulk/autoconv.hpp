#pragma once

#include <algorithm>

namespace libcloudphxx
{
  namespace bulk
  {
    template <typename real_t, class container_t>
    void autoconv(
      container_t drhod_rc_cont,
      container_t drhod_rr_cont,
      const container_t rhod_cont,   
      const container_t rhod_rc_cont,
      const container_t rhod_rr_cont
    )   
    {
      for (auto tup : zip(drhod_rc_cont, drhod_rr_cont, rhod_cont, rhod_rc_cont, rhod_rr_cont))
      {
        real_t
          tmp,
          &drhod_rc = boost::get<0>(tup),
          &drhod_rr = boost::get<1>(tup);
        const real_t
          &rhod     = boost::get<2>(tup),
          &rhod_rc  = boost::get<3>(tup),
          &rhod_rr  = boost::get<4>(tup);

	tmp = std::max( 0., .001 * rhod * (rhod_rc / rhod - .0005)); //should be .001
	drhod_rr += tmp;
	drhod_rc -= tmp;
      }
    }    
  }
};
