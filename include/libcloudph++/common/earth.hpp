#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp>

namespace libcloudphxx
{
  namespace common
  {
    namespace earth
    {
      // acceleration due to gravity
      libcloudphxx_const(si::acceleration, g, 9.81, si::metres_per_second_squared)
    };
  };
};
