#pragma once

#include <algorithm>

namespace libcloudphxx
{
  namespace bulk
  {
    template <class container_t, typename real_t>
    void collect(
      real_t dt,
      const container_t &rhod_cont,   
      container_t &rhod_rc_cont,
      container_t &rhod_rr_cont
    )   
    {
      for (auto tup : zip(rhod_cont, rhod_rc_cont, rhod_rr_cont))
      {
        const real_t
          &rhod = boost::get<0>(tup);
        real_t
          tmp,
          &rhod_rc = boost::get<1>(tup),
          &rhod_rr = boost::get<2>(tup);

        tmp = 2.2 * rhod_rc * ( pow(rhod_rr / rhod, .875));

	rhod_rr += dt * tmp;
	rhod_rc -= dt * tmp;
      }
    }    
  }
};
