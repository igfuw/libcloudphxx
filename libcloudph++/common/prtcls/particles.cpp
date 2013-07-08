#include <iostream>

#include "particles.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      template <int backend, typename real_t>
      particles_proto<real_t> *factory()
      {
        switch (backend)
        {
          case cpp:
	    std::cerr << "allocating CPP..." << std::endl;
	    return new particles<real_t, cpp>();
         
          case omp:
#if defined(_OPENMP)
	    std::cerr << "allocating OpenMP..." << std::endl;
	    return new particles<real_t, omp>();
#else
	    assert(backend != omp && "libcloudph++ was compiled without support for OpenMP");
#endif

          case cuda:
#if defined(CUDA_FOUND) // should be present through CMake's add_definitions()
	    std::cerr << "allocating CUDA..." << std::endl;
	    return new particles<real_t, cuda>();
#else
	    assert(backend != cuda && "libcloudph++ was compiled without support for OpenMP");
#endif

          default:
            assert(false && "unknown backend type");
        }
      }
      
      // explicit instantiation
      template particles_proto<float> *factory<cpp>();
      template particles_proto<float> *factory<omp>();
      template particles_proto<float> *factory<cuda>();

      template particles_proto<double> *factory<cpp>();
      template particles_proto<double> *factory<omp>();
      template particles_proto<double> *factory<cuda>();
    };
  };
};
