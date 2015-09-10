// contains definitions of members of particles_t specialized for multiple GPUs
#include <omp.h>
#include "detail/multiGPU_utils.hpp"

// macro to check for cuda errors, taken from 
// http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define gpuErrchk(ans) { detail::gpuAssert((ans), __FILE__, __LINE__); }

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
      {
         if (code != cudaSuccess) 
         {
            fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
            if (abort) exit(code);
         }
      }
    };

    // constructor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::particles_t(const opts_init_t<real_t> &_opts_init, int dev_id) :
      opts_init(_opts_init)
    {
      int dev_count;
      // TODO: move these sanity checks to sanity_checks?

      // multi_CUDA works only for 2D and 3D
      if(opts_init.nz == 0)
        throw std::runtime_error("multi_CUDA backend works only for 2D and 3D simulations.");

      // get number of available devices
      gpuErrchk(cudaGetDeviceCount(&dev_count)); 
      
      // set number of devices to use
      if(opts_init.dev_count > 0)
      {
        if(dev_count < opts_init.dev_count)
          throw std::runtime_error("number of available GPUs smaller than number of GPUs defined in opts_init");
        else 
          dev_count = opts_init.dev_count;
      }
      // copy dev_count to opts_init for threads to use
      opts_init.dev_count = dev_count;
   
      // check if all GPUs support UVA, TODO: move this to cmake
      for (int i = 0; i < dev_count; ++i)
      {
        // Get device properties
        cudaDeviceProp devProp;
        gpuErrchk(cudaGetDeviceProperties(&devProp, i));
        if(!devProp.unifiedAddressing)
          throw std::runtime_error("One of the GPUs doesn't support Unified Virtual Addressing.");
      }
      
      // resize the pointer vector
      particles.resize(dev_count);

      // assign device to each thread and create particles_t in each
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));

        opts_init_t<real_t> opts_init_tmp(opts_init); // different for each thread
        // modify nx for each device
        opts_init_tmp.nx = detail::get_dev_nx(opts_init, dev_id);

        particles.at(dev_id) = new particles_t<real_t, CUDA>(opts_init_tmp, dev_id);
      }
    }

    // initialisation 
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::init(
      const arrinfo_t<real_t> th,
      const arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3
    )
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->init(th, rv, rhod, courant_1, courant_2, courant_3);
      }
    }

    // time-stepping methods
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::step_sync(
      const opts_t<real_t> &opts_init,
      arrinfo_t<real_t> th,
      arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3,
      const arrinfo_t<real_t> rhod
    )
    {
    }

    template <typename real_t>
    real_t particles_t<real_t, multi_CUDA>::step_async(
      const opts_t<real_t> &opts
    )
    {
    }

    // diagnostic methods
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_sd_conc()
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_mom(const int &k)
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mom(const int &k)
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mass_dens(const real_t&, const real_t&)
    {
    }

    // ...
//</listing>

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_chem(const enum chem_species_t&)
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_rw_ge_rc()
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_RH_ge_Sc()
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_all()
    {
    }

    template <typename real_t>
    real_t* particles_t<real_t, multi_CUDA>::outbuf()
    {
    }

  };
};
