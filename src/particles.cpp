#include <iostream>

#include "../include/lgrngn/particles.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      template <typename real_t>
      particles_proto<real_t> *factory(real_t sd_conc_mean, int nx, int ny, int nz)
      {
#if defined(CUDA_FOUND) // should be present through CMake's add_definitions()
	return new particles<real_t, cuda>(sd_conc_mean, nx, ny, nz);
#elif defined(_OPENMP)
	return new particles<real_t, omp>(sd_conc_mean, nx, ny, nz);
#else
	return new particles<real_t, cpp>(sd_conc_mean, nx, ny, nz);
#endif
      }
      
      // explicit instantiation
      template particles_proto<float>  *factory(float sd_conc_mean, int nx, int ny, int nz);
//      template particles_proto<double> *factory(double sd_conc_mean, int nx, int ny, int nz);
    };
  };
};
