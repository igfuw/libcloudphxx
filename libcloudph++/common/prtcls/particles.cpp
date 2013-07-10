#include <iostream>

#include "particles.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      template <typename real_t>
      particles_proto<real_t> *factory(real_t sd_conc_mean, int nx, int ny, int nz)
      {
        // TODO: provide some controll over the choice

#if defined(CUDA_FOUND) // should be present through CMake's add_definitions()
        std::cerr << "allocating CUDA..." << std::endl;
        return new particles<real_t, cuda>(sd_conc_mean, nx, ny, nz);
#endif
         
#if defined(_OPENMP)
        std::cerr << "allocating OpenMP..." << std::endl;
        return new particles<real_t, omp>(sd_conc_mean, nx, ny, nz);
#endif
 
        std::cerr << "allocating CPP..." << std::endl;
        return new particles<real_t, cpp>(sd_conc_mean, nx, ny, nz);
      }
      
      // explicit instantiation
      template particles_proto<float>  *factory(float  sd_conc_mean, int nx, int ny, int nz);
//      template particles_proto<double> *factory(double sd_conc_mean, int nx, int ny, int nz);
    };
  };
};
