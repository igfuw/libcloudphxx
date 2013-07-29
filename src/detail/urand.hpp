#pragma once

#include "thrust.hpp"

#if defined(__NVCC__)
#  include <curand.h>
#else
#  include <random>
#  include <algorithm>
#endif

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t, int backend>
      class u01
      {
#if !defined(__NVCC__)
	// serial version using C++11's <random>
	using engine_t = std::mt19937;
        using dist_t = std::uniform_real_distribution<real_t>;
	engine_t engine;
	dist_t dist;

	struct fnctr
	{
          engine_t &engine;
          dist_t &dist;
	  real_t operator()() { return dist(engine); }
	};// fnctri({.engine = engine, .dist = dist});

	public:

        // ctor
        u01() : engine(44), dist(0,1) {}

	void generate_n(
	  thrust_device::vector<real_t> &u01, 
	  const thrust_size_t n
	) {
          // note: generate_n copies the third argument!!!
	  std::generate_n(u01.begin(), n, fnctr({.engine = engine, .dist = dist})); 
	}
#endif
      };
 
      template <typename real_t>
      class u01<real_t, cuda>
      {
#if defined(__NVCC__)
	// CUDA parallel version using curand

	// private member fields
	curandGenerator_t gen;
	
	public:

	u01()
	{
	  int status;
	  status = curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MTGP32);
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandCreateGenerator failed"*/);
	  status = curandSetPseudoRandomGeneratorSeed(gen, 44);
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandSetPseudoRandomGeneratorSeed failed"*/);
	}

	~u01()
	{
	  int status = curandDestroyGenerator(gen); 
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandDestroyGenerator failed"*/);
	}

	void generate_n(
	  thrust_device::vector<float> &v, 
	  const thrust_size_t n
	)
	{
	  int status = curandGenerateUniform(gen, thrust::raw_pointer_cast(v.data()), n);
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
	}

	void generate_n(
	  thrust_device::vector<double> &v, 
	  const thrust_size_t n
	)
	{
	  int status = curandGenerateUniformDouble(gen, thrust::raw_pointer_cast(v.data()), n);
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
	}
#endif
      };
    };
  };
};
