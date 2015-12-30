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
      libcloudphxx_const(volume_over_amount,                         R_S_H2O2_K, 13. * 1e-3,       si::cubic_metres / si::moles) 

      libcloudphxx_const(volume_over_amount_over_time, R_S_O3_k0,  2.4 * 1e4 * 1e-3, si::cubic_metres / si::moles / si::seconds) 
      libcloudphxx_const(volume_over_amount_over_time, R_S_O3_k1,  3.7 * 1e5 * 1e-3, si::cubic_metres / si::moles / si::seconds) 
      libcloudphxx_const(volume_over_amount_over_time, R_S_O3_k2,  1.5 * 1e9 * 1e-3, si::cubic_metres / si::moles / si::seconds) 

      //modifications to reaction rates due to temperature, Seinfeld and Pandis table 6.A.7 Sulfur Chemistry
      libcloudphxx_const_derived(si::temperature, dER_H2O2_k,  real_t(-4430.) * si::kelvins)

      libcloudphxx_const_derived(si::temperature, dER_O3_k0, real_t(0.)      * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dER_O3_k1, real_t(-5530.)  * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dER_O3_k2, real_t(-5280.)  * si::kelvins)

      //reaction rate(temperature) for H2O2 reaction
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<volume_square_over_amount_square_over_time, real_t> R_temp_H2O2(
        const quantity<si::temperature, real_t> &T, 
        const quantity<volume_square_over_amount_square_over_time, real_t> &R, 
        const quantity<si::temperature, real_t> &dER
      ) { 
        return (R * exp(dER * (real_t(1.)/T - (real_t(1./298) / si::kelvins))));
      }   
 
      //reaction rate(temperature) for O3 reaction
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<volume_over_amount_over_time, real_t> R_temp_O3(
        const quantity<si::temperature, real_t> &T, 
        const quantity<volume_over_amount_over_time, real_t> &R, 
        const quantity<si::temperature, real_t> &dER
      ) { 
        return (R * exp(dER * (real_t(1.)/T - (real_t(1./298) / si::kelvins))));
      }
   
    };
  };
};
