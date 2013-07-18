#include <iostream>

#include "../include/lgrngn/particles.hpp"

#if defined(_OPENMP)
#  include <omp.h>
#  include "lib.hpp"
#  include <thrust/system/omp/vector.h>
void omp_sanity_check()
{
  if (omp_get_max_threads() == 1) return;
  thrust::omp::vector<int> v(100);
  struct { int operator()(int) { return omp_get_thread_num(); } } thread_id;
  thrust::transform(v.begin(), v.end(), v.begin(), thread_id);
  auto minmax = thrust::minmax_element(v.begin(), v.end());
  assert(*minmax.first != *minmax.second);
}
#endif

namespace libcloudphxx
{
  namespace lgrngn
  {
    // the reasons to have this factory are:
    // - to handle errors like CUDA version not present
    // - to manage 0D/1D/2D/3D parameter defaults
    // - to shorten the code on the caller side
    template <typename real_t>
    particles_proto<real_t> *factory(const int backend, const opts_t<real_t> opts)
    {
      switch (backend)
      {
	case cuda:
#if defined(CUDA_FOUND) // should be present through CMake's add_definitions()
	  return new particles<real_t, cuda>(opts);
#else
	  assert(false && "CUDA backend was not compiled"); throw;
#endif
	case omp:
#if defined(_OPENMP)
          omp_sanity_check();
	  return new particles<real_t, omp>(opts);
#else
	  assert(false && "OpenMP backend was not compiled"); throw;
#endif
	case cpp:
	  return new particles<real_t, cpp>(opts);
	default:
	  assert(false && "unknown backend"); throw;
      }
    }

    template <typename real_t>
    particles_proto<real_t> *factory(
      const int backend, 
      const real_t sd_conc_mean
    )
    {
      opts_t<real_t> opts;
      opts.sd_conc_mean = sd_conc_mean;
      return factory(backend, opts);
    }

    template <typename real_t>
    particles_proto<real_t> *factory(
      const int backend, 
      const real_t sd_conc_mean,
      const int nz, const real_t dz
    )
    {
      opts_t<real_t> opts;
      opts.sd_conc_mean = sd_conc_mean;
      opts.nz = nz;
      opts.dz = dz;
      return factory(backend, opts);
    }

    template <typename real_t>
    particles_proto<real_t> *factory(
      const int backend, 
      const real_t sd_conc_mean,
      const int nx, const int nz,
      const real_t dx, const real_t dz
    )
    {
      opts_t<real_t> opts;
      opts.sd_conc_mean = sd_conc_mean;
      opts.nx = nx;
      opts.nz = nz;
      opts.dx = dx;
      opts.dz = dz;
      return factory(backend, opts);
    }

    template <typename real_t>
    particles_proto<real_t> *factory(
      const int backend, 
      const real_t sd_conc_mean,
      const int nx,
      const int ny,
      const int nz,
      const real_t dx,
      const real_t dy,
      const real_t dz
    )
    {
      opts_t<real_t> opts;
      opts.sd_conc_mean = sd_conc_mean;
      opts.nx = nx;
      opts.ny = ny;
      opts.nz = nz;
      opts.dx = dx;
      opts.dy = dy;
      opts.dz = dz;
      return factory(backend, opts);
    }
    
    // explicit instantiation
    // TODO: what about double precision version?
    template particles_proto<float>  *factory(const int, const float); // 0D
    template particles_proto<float>  *factory(const int, const float, const int, const float); // 1D
    template particles_proto<float>  *factory(const int, const float, const int, const int, const float, const float); // 2D
    template particles_proto<float>  *factory(const int, const float, const int, const int, const int, const float, const float, const float); // 2D
  };
};
