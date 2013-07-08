#include <iostream>

#include "particles.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      template <typename real_t>
      particles_proto<real_t> *factory()
      {
        // TODO: provide some controll over the choice

#if defined(CUDA_FOUND) // should be present through CMake's add_definitions()
        std::cerr << "allocating CUDA..." << std::endl;
        return new particles<real_t, cuda>();
#endif
         
#if defined(_OPENMP)
        std::cerr << "allocating OpenMP..." << std::endl;
        return new particles<real_t, omp>();
#endif
 
        std::cerr << "allocating CPP..." << std::endl;
        return new particles<real_t, cpp>();
      }
      
      // explicit instantiation
      template particles_proto<float> *factory();
      template particles_proto<double> *factory();
    };
  };
};
