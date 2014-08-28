#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp> 
#include <libcloudph++/common/earth.hpp> 

namespace libcloudphxx
{
  namespace common
  {
    typedef divide_typeof_helper<
      si::volume,  
      si::amount
    >::type volume_over_amount;

    typedef divide_typeof_helper<
      volume_over_amount,
      si::time
    >::type volume_over_amount_over_time;

    typedef multiply_typeof_helper<
      typename power_typeof_helper<
        volume_over_amount_over_time,
        static_rational<2>
      >::type,
      si::time
    >::type volume_square_over_amount_square_over_time;

    namespace react
    {
      // reactivity constants for S(IV) -> S(VI) chemical reactions
      // Seinfeld & Pandis 1997
      libcloudphxx_const(volume_square_over_amount_square_over_time, R_S_H2O2_k, 7.5 * 1e7 * 1e-6, si::cubic_metres * si::cubic_metres / si::moles / si::moles / si::seconds) 
      libcloudphxx_const(volume_over_amount, R_S_H2O2_K, 13. * 1e-3, si::cubic_metres / si::moles) 

      libcloudphxx_const(volume_over_amount_over_time, R_S_O3_k0,  2.4 * 1e4 * 1e-3, si::cubic_metres / si::moles / si::seconds) 
      libcloudphxx_const(volume_over_amount_over_time, R_S_O3_k1,  3.7 * 1e5 * 1e-3, si::cubic_metres / si::moles / si::seconds) 
      libcloudphxx_const(volume_over_amount_over_time, R_S_O3_k2,  1.5 * 1e9 * 1e-3, si::cubic_metres / si::moles / si::seconds) 
    };
  };
};
