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
    particles_proto_t<real_t> *factory_helper(const int backend, const opts_init_t<real_t> &opts_init)
    {
      switch (backend)
      {
	case cuda:
#if defined(CUDA_FOUND) // should be present through CMake's add_definitions()
	  return new particles_t<real_t, cuda>(opts_init);
#else
	  assert(false && "CUDA backend was not compiled"); throw; // TODO: convert into exception
#endif
	case omp:
#if defined(_OPENMP)
	  return new particles_t<real_t, omp>(opts_init);
#else
	  assert(false && "OpenMP backend was not compiled"); throw; // TODO: convert into exception
#endif
	case cpp:
	  return new particles_t<real_t, cpp>(opts_init);
	default:
	  assert(false && "unknown backend"); throw; // TODO: convert into exception
      }
    }

    template <typename real_t>
    particles_proto_t<real_t>* factory<real_t>::make(
      const int backend, 
      const opts_init_t<real_t> &opts_init
    ) {
      return factory_helper(backend, opts_init);
    }

    // explicit instantiation
    template struct factory<float>;
    template struct factory<double>;
  };
};
