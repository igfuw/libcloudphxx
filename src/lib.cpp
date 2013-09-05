#include <iostream>

#include <libcloudph++/lgrngn/factory.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // the reasons to have this factory are:
    // - to handle errors like CUDA version not present
    // - to manage 0D/1D/2D/3D parameter defaults
    // - to shorten the code on the caller side
    template <typename real_t>
    particles_proto<real_t> *factory_helper(const int backend, const opts_t<real_t> &opts)
    {
      switch (backend)
      {
	case cuda:
#if defined(CUDA_FOUND) // should be present through CMake's add_definitions()
	  return new particles<real_t, cuda>(opts);
#else
	  assert(false && "CUDA backend was not compiled"); throw; // TODO: convert into exception
#endif
	case omp:
#if defined(_OPENMP)
	  return new particles<real_t, omp>(opts);
#else
	  assert(false && "OpenMP backend was not compiled"); throw; // TODO: convert into exception
#endif
	case cpp:
	  return new particles<real_t, cpp>(opts);
	default:
	  assert(false && "unknown backend"); throw; // TODO: convert into exception
      }
    }

    template <typename real_t>
    particles_proto<real_t>* factory<real_t>::make(
      const int backend, 
      const opts_t<real_t> &opts
    ) {
      return factory_helper(backend, opts);
    }

    // explicit instantiation
    template struct factory<float>;
    template struct factory<double>;
  };
};
