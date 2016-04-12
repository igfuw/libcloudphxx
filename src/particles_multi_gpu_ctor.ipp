// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

// contains definitions of members of particles_t specialized for multiple GPUs
#include <omp.h>
namespace libcloudphxx
{
  namespace lgrngn
  {
    // constructor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::particles_t(opts_init_t<real_t> _opts_init) :
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

      if(glob_opts_init.nx == 0)
        throw std::runtime_error("multi_CUDA doesn't work for 0D setup");

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
      // copy actual dev_count to glob_opts_init
      glob_opts_init.dev_count = dev_count;
   
      // check if all GPUs support UVA
      // TODO: other checks?, see CUDA samples 0_simple/simpleP2P
      for (int i = 0; i < dev_count; ++i)
      {
        // Get device properties
        cudaDeviceProp devProp;
        gpuErrchk(cudaGetDeviceProperties(&devProp, i));
        if(!devProp.unifiedAddressing)
          throw std::runtime_error("All GPUs have to support Unified Virtual Addressing.");
        if(devProp.computeMode != 0)
          throw std::runtime_error("All GPUs have to be in the \"shared\" compute mode.");
      }

      // allow direct memory access between nieghbouring devices
      if(dev_count>1)
      {
        for(int dev_id = 0; dev_id < dev_count; ++dev_id)
        {
          gpuErrchk(cudaSetDevice(dev_id));
          // IDs of devices to the left/right, periodic boundary in x
          const int lft_dev = dev_id > 0 ? dev_id - 1 : glob_opts_init.dev_count - 1,
                    rgt_dev = dev_id < glob_opts_init.dev_count-1 ? dev_id + 1 : 0;

          // if available, allow direct memory access; otherwise copy through host memory will be done
          int can_access_peer;
          gpuErrchk(cudaDeviceCanAccessPeer(&can_access_peer, dev_id, lft_dev));
          if(can_access_peer)
            {gpuErrchk(cudaDeviceEnablePeerAccess(lft_dev, 0));}
          gpuErrchk(cudaDeviceCanAccessPeer(&can_access_peer, dev_id, rgt_dev));
          if(can_access_peer && dev_count > 2)
            {gpuErrchk(cudaDeviceEnablePeerAccess(rgt_dev, 0));}
        }
      }

      #if defined(USE_MPI)
        // initialize mpi with threading support
        const int prov_tlvl = detail::mpi_init(MPI_THREAD_FUNNELED);
        if(prov_tlvl < MPI_THREAD_FUNNELED)
          throw std::runtime_error("MPI was initialized with threading support lower than MPI_THREAD_FUNNELED, multi_CUDA backend won't work");

        // check if it's the main thread of MPI in order to FUNNELED to work
        int main;
        MPI_Is_thread_main(&main);
        if(!main)
          throw std::runtime_error("particles multi_CUDA ctor was called by a thread that is not the main thread of MPI (the mpi_init caller); aborting");
      #endif
      
      // resize the pointer vector
      particles.reserve(dev_count);
      // resize the output buffer
      real_n_cell_tot.resize(n_cell_tot);

      // make opts_init point to global opts init
      this->opts_init = &glob_opts_init;

      // assign device to each thread and create particles_t in each
      int n_x_bfr = 0;
      for(int dev_id = 0; dev_id < dev_count; ++dev_id)
      {
        gpuErrchk(cudaSetDevice(dev_id));
        opts_init_t<real_t> opts_init_tmp(glob_opts_init);

        // adjust opts_init for each device
        if(dev_count > 1)
          n_x_bfr = detail::distmem_opts(opts_init_tmp, dev_id, dev_count); 

        particles.push_back(new particles_t<real_t, CUDA>(opts_init_tmp));
        
        // set n_x_bfr and n_cell_bfr and bcond type for this device 
        particles[dev_id].pimpl->n_x_bfr = n_x_bfr;
        particles[dev_id].pimpl->n_cell_bfr = n_x_bfr * detail::m1(opts_init_tmp.ny) * detail::m1(opts_init_tmp.nz);
        // set distmem types
        if(dev_count > 1)
        {
          if(!particles[dev_id].pimpl->distmem_mpi()) // if there is no MPI copy, set all boundaries to cuda
            particles[dev_id].pimpl->bcond = std::make_pair(detail::distmem_cuda, detail::distmem_cuda);
          else // if there is MPI, set in-node boundaries between devices to cuda
          {
            if(dev_id == 0)
              particles[dev_id].pimpl->bcond.second = detail::distmem_cuda;
            else if(dev_id == dev_count - 1)
              particles[dev_id].pimpl->bcond.first = detail::distmem_cuda;
            else
              particles[dev_id].pimpl->bcond = std::make_pair(detail::distmem_cuda, detail::distmem_cuda);
          }
        }
        // store dev_count in the thread; regular ctor zeroes it
        particles[dev_id].pimpl->opts_init.dev_count = dev_count;
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
      const arrinfo_t<real_t> courant_3,
      const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem
    )
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].init(th, rv, rhod, courant_1, courant_2, courant_3, ambient_chem);

        #pragma omp master
        {
          // exchange domains using mpi; has to be done sequentially here as MPI isn't good with calls from many threads per node
          if(glob_opts_init.dev_count > 1)
          {
            // first node receives first 
            if(particles[0].pimpl->mpi_rank==0)
            {
              particles[0].pimpl->xchng_domains();
              particles[glob_opts_init.dev_count-1].pimpl->xchng_domains();
            }
            else  // other nodes send first
            {
              particles[glob_opts_init.dev_count-1].pimpl->xchng_domains();
              particles[0].pimpl->xchng_domains();
            }
          }
          else
            particles[0].pimpl->xchng_domains();
        }
      }
    }
  };
};
