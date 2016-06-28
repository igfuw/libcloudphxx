#pragma once 

#include "particles.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    // factory will be explicitely instantiated
//<listing>
    template <typename real_t>
    particles_proto_t<real_t> *factory(
      const backend_t, 
      opts_init_t<real_t>
    );
//</listing>
  };
};
