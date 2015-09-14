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
      // compile-time max(1, n), same as in pimpl_ctor...
      int m1(int n) { return n == 0 ? 1 : n; }      
    };

    // constructor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::particles_t(const opts_init_t<real_t> &_opts_init, int dev_id) :
      opts_init(_opts_init),
      n_cell_tot(
        detail::m1(opts_init.nx) *
        detail::m1(opts_init.ny) *
        detail::m1(opts_init.nz)
      ),
      real_n_cell_tot(n_cell_tot)
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
      #pragma omp parallel firstprivate(opts_init) num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));

        // modify nx for each device
        opts_init.nx = detail::get_dev_nx(opts_init, dev_id);

        particles.at(dev_id) = new particles_t<real_t, CUDA>(opts_init, dev_id); // impl stores a copy of opts_init
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
      const opts_t<real_t> &opts,
      arrinfo_t<real_t> th,
      arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3,
      const arrinfo_t<real_t> rhod
    )
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->step_sync(opts, th, rv, courant_1, courant_2, courant_3, rhod);
      }
    }

    template <typename real_t>
    real_t particles_t<real_t, multi_CUDA>::step_async(
      const opts_t<real_t> &opts
    )
    {
      real_t res = 0.;
      #pragma omp parallel reduction(+:res) num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        res = particles[dev_id]->step_async(opts);
      }
      return res;
    }

    // diagnostic methods
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_sd_conc()
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->diag_sd_conc();
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->diag_dry_rng(r_mi, r_mx);
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->diag_wet_rng(r_mi, r_mx);
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_mom(const int &k)
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->diag_dry_mom(k);
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mom(const int &k)
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->diag_wet_mom(k);
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mass_dens(const real_t &a, const real_t &b)
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->diag_wet_mass_dens(a, b);
      }
    }

    // ...
//</listing>

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_chem(const enum chem_species_t &spec)
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->diag_chem(spec);
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_rw_ge_rc()
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->diag_rw_ge_rc();
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_RH_ge_Sc()
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->diag_RH_ge_Sc();
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_all()
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->diag_all();
      }
    }

    template <typename real_t>
    real_t* particles_t<real_t, multi_CUDA>::outbuf()
    {
      #pragma omp parallel num_threads(dev_count)
      {
        const int dev_id = omp_get_thread_num();
        particles[dev_id]->pimpl->fill_outbuf();
      }
      // fill the global output
      int n_cell_bfr;
      for (int dev_id = 0; dev_id < dev_count; ++dev_id)
      {
        n_cell_bfr = dev_id * detail::get_dev_nx(opts_init, 0) * detail::m1(opts_init.ny) * detail::m1(opts_init.nz);
        thrust::copy(
          particles[dev_id]->pimpl->tmp_host_real_cell.begin(),
          particles[dev_id]->pimpl->tmp_host_real_cell.end(),
          thrust::make_permutation_iterator(
            real_n_cell_tot.begin(),
            thrust::make_counting_iterator<int>(n_cell_bfr)
          )
        );
        return &(*(real_n_cell_tot.begin()));
      }
    }
  };
};
