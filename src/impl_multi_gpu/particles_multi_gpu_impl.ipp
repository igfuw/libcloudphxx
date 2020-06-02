// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

// contains definitions of members of particles_t specialized for multiple GPUs
namespace libcloudphxx
{
  namespace lgrngn
  {
    // multi_CUDA pimpl stuff 
    template <typename real_t>
    struct particles_t<real_t, multi_CUDA>::impl
    { 
      std::vector<std::unique_ptr<particles_t<real_t, CUDA> > > particles; // pointer to particles_t on each GPU
      opts_init_t<real_t> glob_opts_init; // global copy of opts_init (threads store their own in impl), 
      const int n_cell_tot;               // total number of cells
      std::vector<real_t> real_n_cell_tot; // vector of the size of the total number of cells to store output

      // cxx threads helper methods
      template<typename F, typename ... Args>
      void mcuda_run(F&& fun, Args&& ... args);

      void step_async_and_copy(
        const opts_t<real_t> &opts,
        const int dev_id,
        std::vector<cudaStream_t> &streams,
        std::vector<cudaEvent_t> &events,
        detail::barrier_t &
      );

      //ctor
      impl(opts_init_t<real_t> _opts_init) :
        glob_opts_init(_opts_init),
        n_cell_tot(
          detail::m1(glob_opts_init.nx) *
          detail::m1(glob_opts_init.ny) *
          detail::m1(glob_opts_init.nz)
        )
      {
        int dev_count;
        // TODO: move these sanity checks to sanity_checks?
        
        if(glob_opts_init.chem_switch) throw std::runtime_error("multi_CUDA is not yet compatible with chemistry. Use other backend or turn off opts_init.chem_switch.");
  
        if(glob_opts_init.nx == 0)
          throw std::runtime_error("multi_CUDA doesn't work for 0D setup.");
  
        if (!(glob_opts_init.x1 > glob_opts_init.x0 && glob_opts_init.x1 <= glob_opts_init.nx * glob_opts_init.dx))
          throw std::runtime_error("!(x1 > x0 & x1 <= min(1,nx)*dx)");
  
        // get number of available devices
        gpuErrchk(cudaGetDeviceCount(&dev_count)); 
        
        // set number of devices to use
        if(glob_opts_init.dev_count > 0)
        {
          if(dev_count < glob_opts_init.dev_count)
            throw std::runtime_error(detail::formatter() << "number of available GPUs (" << dev_count << ") smaller than number of GPUs defined in opts_init ("<< glob_opts_init.dev_count << ")");
          else 
            dev_count = glob_opts_init.dev_count;
        }
  
        if(dev_count > glob_opts_init.nx)
          throw std::runtime_error(detail::formatter() <<"Number of CUDA devices (" << dev_count << ") used is greater than nx (" << glob_opts_init.nx <<")");
  
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
            {
              cudaError_t cuErr = cudaDeviceEnablePeerAccess(lft_dev, 0); 
              if(cuErr != cudaErrorPeerAccessAlreadyEnabled) gpuErrchk(cuErr);
            }
            gpuErrchk(cudaDeviceCanAccessPeer(&can_access_peer, dev_id, rgt_dev));
            if(can_access_peer && dev_count > 2)
            {
              cudaError_t cuErr = cudaDeviceEnablePeerAccess(rgt_dev, 0);
              if(cuErr != cudaErrorPeerAccessAlreadyEnabled) gpuErrchk(cuErr);
            }
          }
        }

        #if defined(USE_MPI)
          // initialize mpi with threading support, TODO: only if it has not been initialize before
          const int prov_tlvl = detail::mpi_init(MPI_THREAD_MULTIPLE);
          if(prov_tlvl < MPI_THREAD_MULTIPLE)
            throw std::runtime_error("MPI was initialized with threading support lower than MPI_THREAD_MULTIPLE, multi_CUDA backend won't work");

          // check if it's the main thread of MPI in order to MULTIPLE to work
          int main;
          MPI_Is_thread_main(&main);
          if(!main)
            throw std::runtime_error("particles multi_CUDA ctor was called by a thread that is not the main thread of MPI (the mpi_init caller); aborting");
        #endif
       
        
        // resize the pointer vector
        particles.reserve(dev_count);
        // resize the output buffer
        real_n_cell_tot.resize(n_cell_tot);
  
        // assign device to each thread and create particles_t in each
        int n_x_bfr = 0;
        for(int dev_id = 0; dev_id < dev_count; ++dev_id)
        {
          gpuErrchk(cudaSetDevice(dev_id));
          opts_init_t<real_t> opts_init_tmp(glob_opts_init);
  
          // adjust opts_init for each device 
          if(dev_count > 1)
            // modify nx for each device
            n_x_bfr = detail::distmem_opts(opts_init_tmp, dev_id, dev_count); 

        //  particles.push_back(new particles_t<real_t, CUDA>(opts_init_tmp, n_x_bfr, glob_opts_init.nx)); // impl stores a copy of opts_init
          particles.emplace_back(std::unique_ptr<particles_t<real_t, CUDA>>(new particles_t<real_t, CUDA>(opts_init_tmp, glob_opts_init.nx))); // impl stores a copy of opts_init

        
          // set n_x_bfr and n_cell_bfr and bcond type for this device 
          particles[dev_id]->pimpl->n_x_bfr = n_x_bfr;
          particles[dev_id]->pimpl->n_cell_bfr = n_x_bfr * detail::m1(opts_init_tmp.ny) * detail::m1(opts_init_tmp.nz);

          // set distmem types: intra-node boundaries between devices to distmem_cuda
          if(dev_count > 1)
          {
            if(dev_id == 0)
              particles[dev_id]->pimpl->bcond.second = detail::distmem_cuda;
            else if(dev_id == dev_count - 1)
              particles[dev_id]->pimpl->bcond.first = detail::distmem_cuda;
            else
              particles[dev_id]->pimpl->bcond = std::make_pair(detail::distmem_cuda, detail::distmem_cuda);

            // if there is no mpi and the boundary is periodic, set outside boundaries to distmem_cuda, since SDs will be copied between different devices on the same node
            // NOTE: we might as well achieve this by replacing all sharedmem boundaries with distmem_cuda?
            if(!particles[dev_id]->pimpl->distmem_mpi() && !opts_init_tmp.open_side_walls)
            {
              if(dev_id == 0)
                particles[dev_id]->pimpl->bcond.first = detail::distmem_cuda;
              else if(dev_id == dev_count - 1)
                particles[dev_id]->pimpl->bcond.second = detail::distmem_cuda;
            }
          }
          // store dev_count in the thread; regular ctor zeroes it
          particles[dev_id]->pimpl->opts_init.dev_count = dev_count;
        }
      }

      // dtor
      ~impl()
      {
        // disable direct memory access between nieghbouring devices (making another instance of prtcls after destruction of first one caused errors due to granting access twice)
        if(glob_opts_init.dev_count>1)
        {
          for(int dev_id = 0; dev_id < glob_opts_init.dev_count; ++dev_id)
          {
            gpuErrchk(cudaSetDevice(dev_id));
            // IDs of devices to the left/right, periodic boundary in x
            const int lft_dev = dev_id > 0 ? dev_id - 1 : glob_opts_init.dev_count - 1,
                      rgt_dev = dev_id < glob_opts_init.dev_count-1 ? dev_id + 1 : 0;

            // if available, allow direct memory access; otherwise copy through host memory will be done
            int can_access_peer;
            gpuErrchk(cudaDeviceCanAccessPeer(&can_access_peer, dev_id, lft_dev));
            if(can_access_peer)
              {gpuErrchk(cudaDeviceDisablePeerAccess(lft_dev));}
            gpuErrchk(cudaDeviceCanAccessPeer(&can_access_peer, dev_id, rgt_dev));
            if(can_access_peer && glob_opts_init.dev_count > 2)
              {gpuErrchk(cudaDeviceDisablePeerAccess(rgt_dev));}
          }
        }
      }
    };

    // run a function concurently on gpus
    template <typename real_t>
    template<typename F, typename ... Args>
    void particles_t<real_t, multi_CUDA>::impl::mcuda_run(F&& fun, Args&& ... args)
    {
      std::vector<std::thread> threads;
      for (int i = 0; i < glob_opts_init.dev_count; ++i)
      {
        threads.emplace_back(
          detail::set_device_and_run, i, 
          std::bind(
            fun,
            particles[i].get(),
            std::forward<Args>(args)...
          )
        );
      }
      for (auto &th : threads) th.join();
    };
  };
};




