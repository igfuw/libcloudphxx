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
        si::amount,  
        si::volume
      >::type,
      si::pressure
    >::type amount_over_volume_over_pressure;

    namespace henry
    {
      // 1e3 (per litre); p_stp (per standard atmosphere)
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_SO2,  real_t(1.23        * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_H2O2, real_t(7.45e4      * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_O3,   real_t(1.13 * 1e-2 * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)

      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_NH3,  real_t(62         * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_HNO3, real_t(2.1 * 1e5  * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_CO2,  real_t(3.4 * 1e-2 * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)
    };
  };
};
