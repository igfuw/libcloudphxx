#include <iostream>

#include <libcloudph++/lgrngn/factory.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // the reasons to have this factory are:
    // - to handle errors like CUDA version not present
    // - to shorten the code on the caller side
    template <typename real_t>
    particles_proto_t<real_t> *factory(const backend_t backend, const opts_init_t<real_t> &opts_init)
    {
      switch (backend)
      {
	case CUDA:
#if defined(CUDA_FOUND) // should be present through CMake's add_definitions()
	  return new particles_t<real_t, CUDA>(opts_init);
#else
	  assert(false && "CUDA backend was not compiled"); throw; // TODO: convert into exception
#endif
	case OpenMP:
#if defined(_OPENMP)
	  return new particles_t<real_t, OpenMP>(opts_init);
#else
	  assert(false && "OpenMP backend was not compiled"); throw; // TODO: convert into exception
#endif
	case serial:
	  return new particles_t<real_t, serial>(opts_init);
	default:
	  assert(false && "unknown backend"); throw; // TODO: convert into exception
      }
    }

    // explicit instantiation
    template particles_proto_t<float> *factory(const backend_t, const opts_init_t<float> &);
    template particles_proto_t<double> *factory(const backend_t, const opts_init_t<double> &);
  };
};
