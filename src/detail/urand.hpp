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
        using engine_t = std::mt19937; // TODO: if real_t = double, use std::mt19937_64
        using dist_u01_t = std::uniform_real_distribution<real_t>;
        using dist_normal01_t = std::normal_distribution<real_t>;
        using dist_size_t = std::uniform_int_distribution<thrust_size_t>;
        engine_t engine;
        dist_u01_t dist_u01;
        dist_normal01_t dist_normal01;
        dist_size_t dist_size;

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
          dist_size_t &dist_size;
          real_t operator()() { return dist_size(engine); }
        };

        public:

        // ctor
        rng(int seed) : engine(seed), dist_u01(0,1), dist_normal01(0,1), dist_size(0, std::numeric_limits<thrust_size_t>::max()) {}

        void reseed(int seed)
        {
          engine.seed(seed);
        }

        void generate_n(
          thrust_device::vector<real_t> &u01, 
          const thrust_size_t n
        ) {
          // note: generate_n copies the third argument!!!
          std::generate_n(u01.begin(), n, fnctr_u01({engine, dist_u01})); // [0,1) range 
        }

        void generate_normal_n(
          thrust_device::vector<real_t> &normal01, 
          const thrust_size_t n
        ) {
          // note: generate_n copies the third argument!!!
          std::generate_n(normal01.begin(), n, fnctr_normal01({engine, dist_normal01})); 
        }

        void generate_n(
          thrust_device::vector<thrust_size_t> &un, 
          const thrust_size_t n
        ) {
          // note: generate_n copies the third argument!!!
          std::generate_n(un.begin(), n, fnctr_un({engine, dist_size})); 
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
          gpuErrchk(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MTGP32));
          gpuErrchk(curandSetPseudoRandomGeneratorSeed(gen, seed));
        }

        void reseed(int seed)
        {
          gpuErrchk(curandSetPseudoRandomGeneratorSeed(gen, seed));
        }

        ~rng()
        {
          gpuErrchk(curandDestroyGenerator(gen));
        }

        void generate_n(
          thrust_device::vector<float> &v, 
          const thrust_size_t n
        )
        {
          gpuErrchk(curandGenerateUniform(gen, thrust::raw_pointer_cast(v.data()), n)); // (0,1] range
          // shift into the expected [0,1) range
          namespace arg = thrust::placeholders;
          thrust::transform(v.begin(), v.begin() + n, v.begin(), float(1) - arg::_1);
        }

        void generate_n(
          thrust_device::vector<double> &v, 
          const thrust_size_t n
        )
        {
          gpuErrchk(curandGenerateUniformDouble(gen, thrust::raw_pointer_cast(v.data()), n)); // (0,1] range
          // shift into the expected [0,1) range
          namespace arg = thrust::placeholders;
          thrust::transform(v.begin(), v.begin() + n, v.begin(), double(1) - arg::_1);
        }

        void generate_normal_n(
          thrust_device::vector<float> &v, 
          const thrust_size_t n
        )
        {
          gpuErrchk(curandGenerateNormal(gen, thrust::raw_pointer_cast(v.data()), n, float(0), float(1)));
        }

        void generate_normal_n(
          thrust_device::vector<double> &v, 
          const thrust_size_t n
        )
        {
          gpuErrchk(curandGenerateNormalDouble(gen, thrust::raw_pointer_cast(v.data()), n, double(0), double(1)));
        }

        void generate_n(
          thrust_device::vector<thrust_size_t> &v, 
          const thrust_size_t n
        )
        {
          gpuErrchk(curandGenerate(gen, thrust::raw_pointer_cast(v.data()), n));
        }
#endif
      };
    };
  };
};
