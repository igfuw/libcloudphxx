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
      impl(const opts_init_t<real_t> &_opts_init) :
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
            throw std::runtime_error(detail::formatter() << "number of available GPUs (" << dev_count << ") smaller than number of GPUs defined in opts_init ("<< glob_opts_init.dev_count << ")");
          else 
            dev_count = glob_opts_init.dev_count;
        }
  
        if(dev_count > glob_opts_init.nx)
          throw std::runtime_error(detail::formatter() <<"Number of CUDA devices (" << dev_count << ") used is greater than nx (" << glob_opts_init.nx <<")");
  
        // copy dev_count to opts_init for threads to use
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
        
        // resize the pointer vector
        particles.reserve(dev_count);
        // resize the output buffer
        real_n_cell_tot.resize(n_cell_tot);
  
        // assign device to each thread and create particles_t in each
        int n_x_bfr;
        for(int dev_id = 0; dev_id < dev_count; ++dev_id)
        {
          gpuErrchk(cudaSetDevice(dev_id));
          opts_init_t<real_t> opts_init_tmp(glob_opts_init);
          n_x_bfr = dev_id * detail::get_dev_nx(glob_opts_init, 0);
  
          if(dev_count > 1)
          {
            // modify nx for each device
            opts_init_tmp.nx = detail::get_dev_nx(glob_opts_init, dev_id);
  
            // adjust x0, x1 for each device
            if(dev_id != 0) opts_init_tmp.x0 = 0.; // TODO: what if x0 greater than domain of first device?
            if(dev_id != dev_count-1) opts_init_tmp.x1 = opts_init_tmp.nx * opts_init_tmp.dx; //TODO: same as above
            else opts_init_tmp.x1 = opts_init_tmp.x1 - n_x_bfr * opts_init_tmp.dx;
  
            // adjust max numer of SDs on each card
            opts_init_tmp.n_sd_max = opts_init_tmp.n_sd_max / dev_count + 1;
          }
        //  particles.push_back(new particles_t<real_t, CUDA>(opts_init_tmp, n_x_bfr, glob_opts_init.nx)); // impl stores a copy of opts_init
          particles.emplace_back(std::unique_ptr<particles_t<real_t, CUDA>>(new particles_t<real_t, CUDA>(opts_init_tmp, n_x_bfr, glob_opts_init.nx))); // impl stores a copy of opts_init
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

    // constructor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::particles_t(const opts_init_t<real_t> &_opts_init) 
    {
      pimpl.reset(new impl(_opts_init));
  
      // make opts_init point to global opts init
      this->opts_init = &(pimpl->glob_opts_init);
    }

    // dtor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::~particles_t() {}

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
      pimpl->mcuda_run(
        &particles_t<real_t, CUDA>::init,
        th, rv, rhod, courant_1, courant_2, courant_3, ambient_chem
      );
    }

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
