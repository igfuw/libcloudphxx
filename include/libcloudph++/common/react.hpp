#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp> 
#include <libcloudph++/common/earth.hpp> 

namespace libcloudphxx
{
  namespace common
  {
    typedef divide_typeof_helper<
      divide_typeof_helper<
        si::volume,  
        si::amount
      >::type,
      si::time
    >::type volume_over_amount_over_time;

    namespace react
    {
      libcloudphxx_const(volume_over_amount_over_time, R_S_O3_k2, 1.5 * 1e9 * 1e-3, si::cubic_metres / si::moles / si::seconds) 
    };
  };
};
