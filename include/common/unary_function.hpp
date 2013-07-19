#pragma once 

#include <cassert>

namespace libcloudphxx
{
  namespace common
  {
    //
    template <typename real_t>
    struct unary_function : std::unary_function<real_t, real_t>
    {
      virtual real_t funval(const real_t) const = 0;
      real_t operator()(const real_t arg) const { return funval(arg); }
 
      // to make Boost ptr_containers work
      virtual unary_function * do_clone() const = 0;
    };

    template <typename real_t>
    inline unary_function<real_t> * new_clone(const unary_function<real_t>& o)
    {
      return o.do_clone();
    }
  };
};
