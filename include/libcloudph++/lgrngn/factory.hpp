#pragma once 

#include "particles.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    // factory will be explicitely instantiated
//<listing>
    template <typename real_t>
    struct factory
    {
      static particles_proto_t<real_t> *make(
        const int backend, 
        const opts_init_t<real_t> &
      );
    };
//</listing>
  };
};
