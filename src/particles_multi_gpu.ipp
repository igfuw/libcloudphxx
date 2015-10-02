// contains definitions of members of particles_t specialized for multiple GPUs
#include <omp.h>
namespace libcloudphxx
{
  namespace lgrngn
  {
    // constructor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::particles_t(const opts_init_t<real_t> &_opts_init, const int &__dev_id, const int &__n_cell_bfr) :
      glob_opts_init(_opts_init),
      n_cell_tot(
        detail::m1(glob_opts_init.nx) *
        detail::m1(glob_opts_init.ny) *
        detail::m1(glob_opts_init.nz)
      )
    {
      int dev_count;
      // TODO: move these sanity checks to sanity_checks?
      
      if(glob_opts_init.src_switch) throw std::runtime_error("multi_CUDA is not yet compatible with source. Use other backend or turn off opts_init.src_switch.");
      if(glob_opts_init.chem_switch) throw std::runtime_error("multi_CUDA is not yet compatible with chemistry. Use other backend or turn off opts_init.chem_switch.");

      // multi_CUDA works only for 2D and 3D
      if(glob_opts_init.nz == 0)
        throw std::runtime_error("multi_CUDA backend works only for 2D and 3D simulations.");

      if (!(glob_opts_init.x1 > glob_opts_init.x0 && glob_opts_init.x1 <= glob_opts_init.nx * glob_opts_init.dx))
        throw std::runtime_error("!(x1 > x0 & x1 <= min(1,nx)*dx)");

      // get number of available devices
      gpuErrchk(cudaGetDeviceCount(&dev_count)); 
      
      // set number of devices to use
      if(glob_opts_init.dev_count > 0)
      {
        if(dev_count < glob_opts_init.dev_count)
          throw std::runtime_error("number of available GPUs smaller than number of GPUs defined in opts_init");
        else 
          dev_count = glob_opts_init.dev_count;
      }
      // copy dev_count to opts_init for threads to use
      glob_opts_init.dev_count = dev_count;
   
      // check if all GPUs support UVA, TODO: move this to cmake
      // TODO: other checks, see CUDA samples
      for (int i = 0; i < dev_count; ++i)
      {
        // Get device properties
        cudaDeviceProp devProp;
        gpuErrchk(cudaGetDeviceProperties(&devProp, i));
        if(!devProp.unifiedAddressing)
          throw std::runtime_error("One of the GPUs doesn't support Unified Virtual Addressing.");
        if(devProp.computeMode != 0)
          throw std::runtime_error("All GPUs used have to be in the \"shared\" compute mode.");
      }
      
      // resize the pointer vector
      particles.reserve(dev_count);
      // resize the output buffer
      real_n_cell_tot.resize(n_cell_tot);

      // make opts_init point to global opts init
      this->opts_init = &glob_opts_init;

      // assign device to each thread and create particles_t in each
      int n_cell_bfr;
      for(int dev_id = 0; dev_id < dev_count; ++dev_id)
      {
        gpuErrchk(cudaSetDevice(dev_id));
        opts_init_t<real_t> opts_init_tmp(glob_opts_init); // firstprivate didn't work
        n_cell_bfr = dev_id * detail::get_dev_nx(glob_opts_init, 0) * detail::m1(glob_opts_init.ny) * detail::m1(glob_opts_init.nz);

        // modify nx for each device
        opts_init_tmp.nx = detail::get_dev_nx(opts_init_tmp, dev_id);
        particles.push_back(new particles_t<real_t, CUDA>(opts_init_tmp, dev_id, n_cell_bfr)); // impl stores a copy of opts_init
      }

      // allow direct memory access between nieghbouring devices
      if(dev_count>1)
      {
        #pragma omp parallel num_threads(dev_count)
        {
          const int dev_id = omp_get_thread_num();
          gpuErrchk(cudaSetDevice(dev_id));
          // to the left
          if(dev_id != 0)
            {gpuErrchk(cudaDeviceEnablePeerAccess(dev_id-1, 0));}
          else
            {gpuErrchk(cudaDeviceEnablePeerAccess(dev_count-1, 0));}
          // to the right
          if(dev_count > 2)
          {
            if(dev_id != dev_count-1)
              {gpuErrchk(cudaDeviceEnablePeerAccess(dev_id+1, 0));} 
            else
              {gpuErrchk(cudaDeviceEnablePeerAccess(0, 0));}
          }
        }
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
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].init(th, rv, rhod, courant_1, courant_2, courant_3);
      }
    }
  };
};
