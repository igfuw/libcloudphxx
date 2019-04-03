#pragma once 

#include <cassert>
#include <boost/config.hpp>

namespace libcloudphxx
{
  namespace common
  {
    //
    template <typename real_t>
    struct unary_function : std::unary_function<real_t, real_t>
    {
      BOOST_GPU_ENABLED
      virtual real_t funval(const real_t) const = 0;
      BOOST_GPU_ENABLED
      real_t operator()(const real_t arg) const { return funval(arg); }
    };
  };
};
