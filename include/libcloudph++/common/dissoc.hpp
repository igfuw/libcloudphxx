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
      // TODO temperature dependence (provided values are for T=298K)

      // bec. there are so few water ions instead of K we have K [H2O]
      // see Seinfeld & Pandis p 345
      libcloudphxx_const(amount_square_over_volume_square, K_H2O, 1e-14 * 1e6 , si::moles * si::moles / si::cubic_metres / si::cubic_metres)

      libcloudphxx_const(amount_over_volume, K_SO2,  1.3e-2 * 1e3, si::moles / si::cubic_metres)
      libcloudphxx_const(amount_over_volume, K_HSO3, 6.6e-8 * 1e3, si::moles / si::cubic_metres)

      // assumes no not-dissociated H2SO4 (see Seinfeld & Pandiss p. 388)
      libcloudphxx_const(amount_over_volume, K_HSO4, 1.2e-2 * 1e3, si::moles / si::cubic_metres) 
    };
  };
};
