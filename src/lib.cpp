#include <iostream>
#include <exception>
#include <libcloudph++/lgrngn/factory.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // the reasons to have this factory are:
    // - to handle errors like CUDA version not present
    // - to shorten the code on the caller side
    template <typename real_t>
    particles_proto_t<real_t> *factory(const backend_t backend, opts_init_t<real_t> opts_init)
    {
      if(backend != multi_CUDA) opts_init.dev_count = 0; // override user-defined dev_count if not using multi_CUDA

      switch (backend)
      {
	case multi_CUDA:
#if defined(CUDA_FOUND) // should be present through CMake's add_definitions(), TODO: some other check in CMake?
	  return new particles_t<real_t, multi_CUDA>(opts_init);
#else
          throw std::runtime_error("multi_CUDA backend was not compiled");
#endif
	case CUDA:
#if defined(CUDA_FOUND) // should be present through CMake's add_definitions()
	  return new particles_t<real_t, CUDA>(opts_init);
#else
          throw std::runtime_error("CUDA backend was not compiled");
#endif
	case OpenMP:
#if defined(_OPENMP)
	  return new particles_t<real_t, OpenMP>(opts_init);
#else
          throw std::runtime_error("OpenMP backend was not compiled"); 
#endif
	case serial:
	  return new particles_t<real_t, serial>(opts_init);
	default:
          throw std::runtime_error("unknown backend"); 
      }
    }

    // explicit instantiation
    template particles_proto_t<float> *factory(const backend_t, opts_init_t<float>);
    template particles_proto_t<double> *factory(const backend_t, opts_init_t<double>);
  };
};
