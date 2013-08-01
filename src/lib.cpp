#include <iostream>

#include <libcloudph++/lgrngn/particles.hpp>

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
/*    
      const real_t sd_conc_mean,
      typename opts_t<real_t>::dry_distros_t dry_distros,
      const int n1, const real_t d1,
      const int n2, const real_t d2,
      const int n3, const real_t d3
    ) {
      opts_t<real_t> opts;
      opts.sd_conc_mean = sd_conc_mean;
      opts.dry_distros = std::move(dry_distros);

      if (n1 == -1 && n2 == -1 && n3 == -1) // 0D
      {
	opts.nz = 0;  opts.dz = 1;
	opts.nx = 0;  opts.dx = 1; 
	opts.ny = 0;  opts.dy = 1; 
      }
      else if (n1 > 0 && n2 == -1 && n3 == -1) // 1D
      {
	opts.nz = n1; opts.dz = d1;
	opts.nx = 0;  opts.dx = 1; 
	opts.ny = 0;  opts.dy = 1; 
      }
      else if (n1 > 0 && n2 > 0 && n3 == -1) // 2D
      {
	opts.nx = n1; opts.dx = d1;
	opts.nz = n2; opts.dz = d2;
	opts.ny = 0;  opts.dy =  1;
      }
      else if (n1 > 0 && n2 > 0 && n3 > 0) // 3D
      {
	opts.nx = n1; opts.dx = d1;
	opts.ny = n2; opts.dy = d2;
	opts.nz = n3; opts.dz = d3;
      }
      else assert(false);
      return factory_helper(backend, opts);
    }
*/

    // explicit instantiation
    template struct factory<float>;
    template struct factory<double>;
  };
};
