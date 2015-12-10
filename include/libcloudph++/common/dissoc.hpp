#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp>

namespace libcloudphxx
{
  namespace common
  {
    typedef divide_typeof_helper<si::amount, si::volume>::type amount_over_volume;
    typedef power_typeof_helper<amount_over_volume, static_rational<2> >::type amount_square_over_volume_square;

    namespace dissoc
    {
      // dissociation constants for chemical compounds of droplets

      // pure water 
      // bec. there are so few water ions instead of K we have K [H2O]
      // see Seinfeld & Pandis p 345
      libcloudphxx_const(amount_square_over_volume_square, K_H2O, 1e-14 * 1e6 , si::moles * si::moles / si::cubic_metres / si::cubic_metres)

      //SO2 * H20 -> HSO3, SO3
      libcloudphxx_const(amount_over_volume, K_SO2,  1.3e-2 * 1e3, si::moles / si::cubic_metres)
      libcloudphxx_const(amount_over_volume, K_HSO3, 6.6e-8 * 1e3, si::moles / si::cubic_metres)

      //H2SO4 -> HSO4, SO4
      // assumes no not-dissociated H2SO4 (see Seinfeld & Pandiss p. 388)
      libcloudphxx_const(amount_over_volume, K_HSO4, 1.2e-2 * 1e3, si::moles / si::cubic_metres) 

      //CO2 * H20 -> HCO3, CO3
      libcloudphxx_const(amount_over_volume, K_CO2,  4.3e-7   * 1e3, si::moles / si::cubic_metres)
      libcloudphxx_const(amount_over_volume, K_HCO3, 4.68e-11 * 1e3, si::moles / si::cubic_metres)

      //NH3 * H20 -> NH4
      libcloudphxx_const(amount_over_volume, K_NH3,  1.7e-5 * 1e3, si::moles / si::cubic_metres)

      //HNO3 -> NO3
      libcloudphxx_const(amount_over_volume, K_HNO3,  15.4 * 1e3, si::moles / si::cubic_metres)

      //modifications to dissociation constants due to temperature, Table 4 in Kreidenweiss et al 2003
      libcloudphxx_const_derived(si::temperature, dKR_CO2,  real_t(-1000.) * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dKR_HCO3, real_t(-1760.) * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dKR_SO2,  real_t(1960.)  * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dKR_HSO3, real_t(1500.)  * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dKR_NH3,  real_t(-450.)  * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dKR_HNO3, real_t(8700.)  * si::kelvins)
      // Seinfeld and Pandis table 6.A.1 Equilibrum Reactions
      libcloudphxx_const_derived(si::temperature, dKR_HSO4, real_t(2720.)  * si::kelvins)

      //dissociation constant(temperature)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<amount_over_volume, real_t> K_temp(
        const quantity<si::temperature, real_t> &T,
        const quantity<amount_over_volume, real_t> &K,
        const quantity<si::temperature, real_t> &dKR
      ) {
        return (K * exp(dKR * (real_t(1.)/T - (real_t(1./298) / si::kelvins))));
      }

    };
  };
};
