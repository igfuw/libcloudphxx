#include <iostream>

#include "../include/lgrngn/particles.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      template <typename real_t>
      particles_proto<real_t> *factory(int backend, real_t sd_conc_mean, int nx, int ny, int nz)
      {
        switch (backend)
        {
          case cuda:
#if defined(CUDA_FOUND) // should be present through CMake's add_definitions()
            return new particles<real_t, cuda>(sd_conc_mean, nx, ny, nz);
#else
            assert(false && "CUDA backend was not compiled"); throw;
#endif
          case omp:
#if defined(_OPENMP)
            return new particles<real_t, omp>(sd_conc_mean, nx, ny, nz);
#else
            assert(false && "OpenMP backend was not compiled"); throw;
#endif
          case cpp:
	    return new particles<real_t, cpp>(sd_conc_mean, nx, ny, nz);
          default:
            assert(false && "unknown backend"); throw;
        }
      }
      
      // explicit instantiation
      template particles_proto<float>  *factory(int backend, float sd_conc_mean, int nx, int ny, int nz);
//      template particles_proto<double> *factory(int backend, double sd_conc_mean, int nx, int ny, int nz);
    };
  };
};
