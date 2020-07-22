#pragma once

#include "thrust.hpp"

#if defined(__NVCC__)
#  include <curand.h>
#  include <limits>
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
      class rng
      {
#if !defined(__NVCC__)
	// serial version using C++11's <random>
	using engine_t = std::mt19937;
        using dist_u01_t = std::uniform_real_distribution<real_t>;
        using dist_normal01_t = std::normal_distribution<real_t>;
        using dist_un_t = std::uniform_int_distribution<unsigned int>;
	engine_t engine;
	dist_u01_t dist_u01;
	dist_normal01_t dist_normal01;
	dist_un_t dist_un;

	struct fnctr_u01
	{
          engine_t &engine;
          dist_u01_t &dist_u01;
	  real_t operator()() { return dist_u01(engine); }
	};

	struct fnctr_normal01
	{
          engine_t &engine;
          dist_normal01_t &dist_normal01;
	  real_t operator()() { return dist_normal01(engine); }
	};

	struct fnctr_un
	{
          engine_t &engine;
          dist_un_t &dist_un;
	  real_t operator()() { return dist_un(engine); }
	};

	public:

        // ctor
        rng(int seed) : engine(seed), dist_u01(0,1), dist_normal01(0,1), dist_un(0, std::numeric_limits<unsigned int>::max()) {}

	void generate_n(
	  thrust_device::vector<real_t> &u01, 
	  const thrust_size_t n
	) {
          // note: generate_n copies the third argument!!!
	  std::generate_n(u01.begin(), n, fnctr_u01({engine, dist_u01})); 
	}

	void generate_normal_n(
	  thrust_device::vector<real_t> &normal01, 
	  const thrust_size_t n
	) {
          // note: generate_n copies the third argument!!!
	  std::generate_n(normal01.begin(), n, fnctr_normal01({engine, dist_normal01})); 
	}

	void generate_n(
	  thrust_device::vector<unsigned int> &un, 
	  const thrust_size_t n
	) {
          // note: generate_n copies the third argument!!!
	  std::generate_n(un.begin(), n, fnctr_un({engine, dist_un})); 
	}
#endif
      };
 
      template <typename real_t>
      class rng<real_t, CUDA>
      {
#if defined(__NVCC__)
	// CUDA parallel version using curand

	// private member fields
	curandGenerator_t gen;
	
	public:

	rng(int seed)
	{
          {
	    int status = curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MTGP32);
	    assert(status == CURAND_STATUS_SUCCESS /* && "curandCreateGenerator failed"*/);
	    _unused(status);
          }
          {
	    int status = curandSetPseudoRandomGeneratorSeed(gen, seed);
	    assert(status == CURAND_STATUS_SUCCESS /* && "curandSetPseudoRandomGeneratorSeed failed"*/);
            _unused(status);
	  }
        }

	~rng()
	{
	  int status = curandDestroyGenerator(gen); 
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandDestroyGenerator failed"*/);
          _unused(status);
	}

	void generate_n(
	  thrust_device::vector<float> &v, 
	  const thrust_size_t n
	)
	{
	  int status = curandGenerateUniform(gen, thrust::raw_pointer_cast(v.data()), n);
          assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
          _unused(status);

	}

	void generate_n(
	  thrust_device::vector<double> &v, 
	  const thrust_size_t n
	)
	{
	  int status = curandGenerateUniformDouble(gen, thrust::raw_pointer_cast(v.data()), n);
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
          _unused(status);
	}

	void generate_normal_n(
	  thrust_device::vector<float> &v, 
	  const thrust_size_t n
	)
	{
	  int status = curandGenerateNormal(gen, thrust::raw_pointer_cast(v.data()), n, float(0), float(1));
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
          _unused(status);
	}

	void generate_normal_n(
	  thrust_device::vector<double> &v, 
	  const thrust_size_t n
	)
	{
	  int status = curandGenerateNormalDouble(gen, thrust::raw_pointer_cast(v.data()), n, double(0), double(1));
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
          _unused(status);
	}

	void generate_n(
	  thrust_device::vector<unsigned int> &v, 
	  const thrust_size_t n
	)
	{
	  int status = curandGenerate(gen, thrust::raw_pointer_cast(v.data()), n);
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
          _unused(status);
	}
#endif
      };
    };
  };
};
