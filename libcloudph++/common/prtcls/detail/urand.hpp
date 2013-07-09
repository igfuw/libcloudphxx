#pragma once

#include "thrust.hpp"

#if (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA)
#else
#  include <thrust/random.h>
#endif

namespace libcloudphxx
{
namespace common
{
namespace prtcls
{
namespace detail
{

template <typename real_t>
void urand(
  thrust::device_vector<real_t> &u01, 
  const thrust_size_t n,
  real_t *seed
)
{
#if (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CPP)
  struct rng
  {
    // member fields
    real_t *seed;
    thrust::random::taus88 engine;
    thrust::uniform_real_distribution<real_t> dist;

    // ctor
    rng(real_t *seed) : seed(seed), engine(*seed), dist(0, 1) {}

    // dtor
    ~rng() { *seed = this->operator()(); }

    // overloaded op invoked by generate()
    real_t operator()() { return dist(engine); }
  };

  thrust::generate_n(u01.begin(), n, rng(seed));
#endif
}

};
};
};
};
