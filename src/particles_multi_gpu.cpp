// contains definitions of members of particles_t specialized for multiple GPUs
#include <omp.h>
#include "detail/multiGPU_utils.hpp"

// macro to check for cuda errors, taken from 
// move it to utils...
// http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define gpuErrchk(ans) { detail::gpuAssert((ans), __FILE__, __LINE__); }

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
      {
         if (code != cudaSuccess) 
         {
            fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
            if (abort) exit(code);
         }
      }

      // max(1, n)
      int m1(int n) { return n == 0 ? 1 : n; }
    };

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
//printf("dev count %d\n", dev_count);
   
      // check if all GPUs support UVA, TODO: move this to cmake
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
//      particles.resize(dev_count);
      particles.reserve(dev_count);
      // resize the output buffer
      real_n_cell_tot.resize(n_cell_tot);

      // make opts_init point to global opts init
      this->opts_init = &glob_opts_init;

      // assign device to each thread and create particles_t in each
//      #pragma omp parallel num_threads(dev_count)
      int n_cell_bfr;
      for(int dev_id = 0; dev_id < dev_count; ++dev_id)
      {
//printf("get thread num\n");
//        const int dev_id = omp_get_thread_num();
//printf("set dev to %d\n", dev_id);
        gpuErrchk(cudaSetDevice(dev_id));

//printf("init temp opts_init\n");
        opts_init_t<real_t> opts_init_tmp(glob_opts_init); // firstprivate didn't work

//printf("calc n_cell_bfr\n");
        /*const int */n_cell_bfr = dev_id * detail::get_dev_nx(glob_opts_init, 0) * detail::m1(glob_opts_init.ny) * detail::m1(glob_opts_init.nz);

//printf("n cell bfr %d\n", n_cell_bfr);
        // modify nx for each device
        opts_init_tmp.nx = detail::get_dev_nx(opts_init_tmp, dev_id);
//printf("local nx %d\n", opts_init_tmp.nx);

        particles.push_back(new particles_t<real_t, CUDA>(opts_init_tmp, dev_id, n_cell_bfr)); // impl stores a copy of opts_init
//        particles.at(dev_id) = boost::make_shared<particles_t<real_t, CUDA> >(opts_init_tmp, dev_id, n_cell_bfr); // impl stores a copy of opts_init
      }
      // allow direct memory access between nieghbouring devices
      // and create stream for each device
      if(dev_count>1)
      {
        streams = new cudaStream_t[dev_count+1]; // stream[0] left for host - see CUDA example 0_simple/UnifiedMemoryStreams
        events = new cudaEvent_t[dev_count+1]; // stream[0] left for host - see CUDA example 0_simple/UnifiedMemoryStreams
        cudaStreamCreate(&streams[0]);
        cudaEventCreateWithFlags (&events[0], cudaEventDisableTiming );
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
          // create a stream
          gpuErrchk(cudaStreamCreate(&streams[dev_id+1]));
          // create an event
          cudaEventCreateWithFlags(&events[dev_id+1], cudaEventDisableTiming );
        }
      }
    }
    // TODO: move methods to a separate file

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
        #pragma omp critical
        {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        printf("call init on dev %d\n", dev_id);
        particles[dev_id].init(th, rv, rhod, courant_1, courant_2, courant_3);
        }
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
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        #pragma omp critical
        {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        printf("call sync on dev %d\n", dev_id);
        particles[dev_id].step_sync(opts, th, rv, courant_1, courant_2, courant_3, rhod);
        }
      }
    }

    template <typename real_t>
    real_t particles_t<real_t, multi_CUDA>::step_async(
      const opts_t<real_t> &opts
    )
    {
      real_t res = 0.;
      #pragma omp parallel reduction(+:res) num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));

        #pragma omp critical
        {
        // do step async on each device
        printf("call async on dev %d\n", dev_id);
        res = particles[dev_id].step_async(opts);
        }

        // --- copy advected SDs to other devices ---
        if(opts.adve && glob_opts_init.dev_count>1)
        {
          printf("start of copying on dev %d\n", dev_id);
          namespace arg = thrust::placeholders;
          typedef unsigned long long n_t; // TODO: same typedef is in impl struct !! TODO: particles::impl::n_t ? 

          // helper aliases
          const unsigned int &lft_count(particles[dev_id].pimpl->lft_count);
          const unsigned int &rgt_count(particles[dev_id].pimpl->rgt_count);
          thrust_size_t &n_part(particles[dev_id].pimpl->n_part);
          thrust_size_t &n_part_old(particles[dev_id].pimpl->n_part_old);
          thrust_device::vector<real_t> &x(particles[dev_id].pimpl->x);
          thrust_device::vector<real_t> &y(particles[dev_id].pimpl->y);
          thrust_device::vector<real_t> &z(particles[dev_id].pimpl->z);
          thrust_device::vector<real_t> &rd3(particles[dev_id].pimpl->rd3);
          thrust_device::vector<real_t> &rw2(particles[dev_id].pimpl->rw2);
          thrust_device::vector<real_t> &kpa(particles[dev_id].pimpl->kpa);
          thrust_device::vector<real_t> &out_real_bfr(particles[dev_id].pimpl->out_real_bfr);
          thrust_device::vector<real_t> &in_real_bfr(particles[dev_id].pimpl->in_real_bfr);
          thrust_device::vector<n_t> &n(particles[dev_id].pimpl->n);
          thrust_device::vector<n_t> &out_n_bfr(particles[dev_id].pimpl->out_n_bfr);
          thrust_device::vector<n_t> &in_n_bfr(particles[dev_id].pimpl->in_n_bfr);
          // i and k must have not changed since impl->bcnd !!
          const thrust_device::vector<thrust_size_t> &lft_id(particles[dev_id].pimpl->i);
          const thrust_device::vector<thrust_size_t> &rgt_id(particles[dev_id].pimpl->k);

          const int lft_dev = dev_id > 0 ? dev_id - 1 : glob_opts_init.dev_count - 1, // periodic boundary in x
                    rgt_dev = dev_id < glob_opts_init.dev_count-1 ? dev_id + 1 : 0; // periodic boundary in x

          // prepare buffer with n_t to be copied left
          // TODO: use serialization to pack the buffers (e.g. from boost)?
          thrust::copy(
            thrust::make_permutation_iterator(n.begin(), lft_id.begin()),
            thrust::make_permutation_iterator(n.begin(), lft_id.begin()) + lft_count,
            out_n_bfr.begin()
          );

          #pragma omp critical
          {
            printf("przed copy dev %d lft_count %d rgt count %d npart %d lft_dev %d rgt_dev %d n_part %d\n", dev_id, lft_count, rgt_count, n_part, lft_dev, rgt_dev, int(n_part));
            printf("lft_id\n");
            debug::print(lft_id.begin(), lft_id.begin() + lft_count);
            printf("rgt_id\n");
            debug::print(rgt_id.begin(), rgt_id.begin() + rgt_count);
            printf("n\n");
            debug::print(n.begin(), n.end());
            printf("x\n");
            debug::print(x.begin(), x.end());
  
            printf("out_n_bfr left\n");
            debug::print(out_n_bfr.begin(), out_n_bfr.begin() + lft_count);
          }

          printf("copy 2 on dev %d\n", dev_id);
          // start async copy of n buffer to the left
          gpuErrchk(cudaMemcpyPeerAsync(
            particles[lft_dev].pimpl->in_n_bfr.data().get(), lft_dev,  //dst
            out_n_bfr.data().get(), dev_id,                             //src 
            lft_count * sizeof(n_t),                                    //no of bytes
            streams[dev_id+1]                                             //best performance if stream belongs to src
          ));
          // record beginning of copying
          gpuErrchk(cudaEventRecord(events[dev_id+1], streams[dev_id+1]));
          // barrier to make sure that all devices started copying
          #pragma omp barrier

          printf("copy 3 on dev %d\n", dev_id);
          // adjust x of prtcls to be sent left to match new device's domain
          thrust::transform(
            thrust::make_permutation_iterator(x.begin(), lft_id.begin()),
            thrust::make_permutation_iterator(x.begin(), lft_id.begin()) + lft_count,
            thrust::make_permutation_iterator(x.begin(), lft_id.begin()), // in place
            arg::_1 + particles[lft_dev].opts_init->x1 - particles[dev_id].opts_init->x0   // operation
          );

          printf("copy 4 on dev %d\n", dev_id);
          // prepare the real_t buffer for copy left
          thrust_device::vector<real_t> * real_t_vctrs[] = {&rd3, &rw2, &kpa, &x, &z, &y};
          const int real_vctrs_count = glob_opts_init.ny == 0 ? 5 : 6;
          for(int i = 0; i < real_vctrs_count; ++i)
            thrust::copy(
              thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), lft_id.begin()),
              thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), lft_id.begin()) + lft_count,
              out_real_bfr.begin() + i * lft_count
            );

          #pragma omp critical
          {
            printf("po zmianie lewych x dev %d lft_count %d rgt count %d npart %d lft_dev %d rgt_dev %d n_part %d\n", dev_id, lft_count, rgt_count, n_part, lft_dev, rgt_dev, int(n_part));
            printf("x\n");
            debug::print(x.begin(), x.end());
            printf("out real bfr na lewo\n");
            debug::print(out_real_bfr.begin(), out_real_bfr.begin() + real_vctrs_count * lft_count);
          }
          

          printf("copy 5 on dev %d\n", dev_id);
          // wait for the copy of n from right into current device to finish
          gpuErrchk(cudaStreamWaitEvent(streams[dev_id+1], events[rgt_dev+1], 0));
          printf("copy 6 on dev %d\n", dev_id);
          // unpack the n buffer sent to this device from right
          int n_copied = particles[rgt_dev].pimpl->lft_count;
          n_part_old = n_part;
          n_part += n_copied;
          n.resize(n_part);
          thrust::copy(in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);
          #pragma omp critical
          {
            printf("po przyjeciu n z prawej dev %d lft_count %d rgt count %d npart %d lft_dev %d rgt_dev %d n_part %d\n", dev_id, lft_count, rgt_count, n_part, lft_dev, rgt_dev, int(n_part));
            printf("n copied from right %d\n", n_copied);
            printf("in_n_bfr z prawej\n");
            debug::print(in_n_bfr.begin(), in_n_bfr.begin() + n_copied);
            printf("n\n");
            debug::print(n.begin(), n.end());
          }

          printf("copy 7 on dev %d\n", dev_id);
          // start async copy of real buffer to the left; same stream as n_bfr - will start only if previous copy finished
          cudaMemcpyPeerAsync(
            particles[lft_dev].pimpl->in_real_bfr.data().get(), lft_dev,  //dst
            out_real_bfr.data().get(), dev_id,                             //src 
            real_vctrs_count * lft_count * sizeof(real_t),                 //no of bytes
            streams[dev_id+1]                                                //best performance if stream belongs to src
          );
          // record beginning of copying
          gpuErrchk(cudaEventRecord(events[dev_id+1], streams[dev_id+1]));
          // barrier to make sure that all devices started copying
          #pragma omp barrier

          printf("copy 8 on dev %d\n", dev_id);
          // prepare buffer with n_t to be copied right
          thrust::copy(
            thrust::make_permutation_iterator(n.begin(), rgt_id.begin()),
            thrust::make_permutation_iterator(n.begin(), rgt_id.begin()) + rgt_count,
            out_n_bfr.begin()
          );
          #pragma omp critical
          {
            printf("gotowy bufor n prawo dev %d lft_count %d rgt count %d npart %d lft_dev %d rgt_dev %d n_part %d\n", dev_id, lft_count, rgt_count, n_part, lft_dev, rgt_dev, int(n_part));
  
            printf("out_n_bfr na prawo\n");
            debug::print(out_n_bfr.begin(), out_n_bfr.begin() + rgt_count);
          }

          printf("copy 9 on dev %d\n", dev_id);
          // adjust x of prtcls to be sent right to match new device's domain
          thrust::transform(
            thrust::make_permutation_iterator(x.begin(), rgt_id.begin()),
            thrust::make_permutation_iterator(x.begin(), rgt_id.begin()) + rgt_count,
            thrust::make_permutation_iterator(x.begin(), rgt_id.begin()), // in place
            arg::_1 + particles[rgt_dev].opts_init->x0 - particles[dev_id].opts_init->x1   // operation
          );

          // wait for the copy of real from right into current device to finish
          gpuErrchk(cudaStreamWaitEvent(streams[dev_id+1], events[rgt_dev+1], 0));

          // unpack the real buffer sent to this device from right
          for(int i = 0; i < real_vctrs_count; ++i)
          {
            real_t_vctrs[i]->resize(n_part);
            thrust::copy(in_real_bfr.begin() + i * n_copied, in_real_bfr.begin() + (i+1) * n_copied, real_t_vctrs[i]->begin() + n_part_old);
          }
          #pragma omp critical
          {
            printf("po zmianie prawych x i rozpakowaniu bufora z prawej dev %d lft_count %d rgt count %d npart %d lft_dev %d rgt_dev %d n_part %d\n", dev_id, lft_count, rgt_count, n_part, lft_dev, rgt_dev, int(n_part));
            printf("x\n");
            debug::print(x.begin(), x.end());
            printf("in real bfr z prawej\n");
            debug::print(in_real_bfr.begin(), in_real_bfr.begin() + real_vctrs_count * n_copied);
          }

          printf("copy 10 on dev %d\n", dev_id);
          // start async copy of n buffer to the right
          cudaMemcpyPeerAsync(
            particles[rgt_dev].pimpl->in_n_bfr.data().get(), rgt_dev,  //dst
            out_n_bfr.data().get(), dev_id,                             //src 
            rgt_count * sizeof(n_t),                                    //no of bytes
            streams[dev_id+1]                                             //best performance if stream belongs to src
          );
          // record beginning of copying
          gpuErrchk(cudaEventRecord(events[dev_id+1], streams[dev_id+1]));
          // barrier to make sure that all devices started copying
          #pragma omp barrier

          // prepare the real_t buffer for copy to the right
          for(int i = 0; i < real_vctrs_count; ++i)
            thrust::copy(
              thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), rgt_id.begin()),
              thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), rgt_id.begin()) + rgt_count,
              out_real_bfr.begin() + i * rgt_count
            );

          // wait for the copy of n from left into current device to finish
          gpuErrchk(cudaStreamWaitEvent(streams[dev_id+1], events[lft_dev+1], 0));
          // unpack the n buffer sent to this device from left
          n_copied = particles[lft_dev].pimpl->rgt_count;
          printf("n copied from left %d\n", n_copied);
          n_part_old = n_part;
          n_part += n_copied;
          n.resize(n_part);
          thrust::copy(in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);

          printf("copy 11 on dev %d\n", dev_id);
          // start async copy of real buffer to the right
          cudaMemcpyPeerAsync(
            particles[rgt_dev].pimpl->in_real_bfr.data().get(), rgt_dev,  //dst
            out_real_bfr.data().get(), dev_id,                             //src 
            real_vctrs_count * rgt_count * sizeof(real_t),                 //no of bytes
            streams[dev_id+1]                                                //best performance if stream belongs to src
          );
          // record beginning of copying
          gpuErrchk(cudaEventRecord(events[dev_id+1], streams[dev_id+1]));
          // barrier to make sure that all devices started copying
          #pragma omp barrier

          printf("copy 12 on dev %d\n", dev_id);
          // flag SDs sent left/right for removal
          thrust::copy(
            thrust::make_constant_iterator<n_t>(0),
            thrust::make_constant_iterator<n_t>(0) + lft_count,
            thrust::make_permutation_iterator(n.begin(), lft_id.begin())
          );
          thrust::copy(
            thrust::make_constant_iterator<n_t>(0),
            thrust::make_constant_iterator<n_t>(0) + rgt_count,
            thrust::make_permutation_iterator(n.begin(), rgt_id.begin())
          );
          
          // wait for the copy of real from left into current device to finish
          gpuErrchk(cudaStreamWaitEvent(streams[dev_id+1], events[lft_dev+1], 0));
          printf("copy 13 on dev %d\n", dev_id);

          // unpack the real buffer sent to this device from left
          for(int i = 0; i < real_vctrs_count; ++i)
          {
            real_t_vctrs[i]->resize(n_part);
            thrust::copy(in_real_bfr.begin() + i * n_copied, in_real_bfr.begin() + (i+1) * n_copied, real_t_vctrs[i]->begin() + n_part_old);
          }

          // particles are not sorted now
          particles[dev_id].pimpl->sorted = false;          
          printf("copy end on dev %d\n", dev_id);
        }

        // finalize async, same as in impl_step - TODO: move to a single function...
        if(glob_opts_init.dev_count>1)
        {   
          // recycling out-of-domain/invalidated particles 
          // currently DISABLED
          thrust_size_t n_rcyc = 0;//pimpl->rcyc();
          // TODO: ! if we do not recycle, we should remove them to care for out-od-domain advection after sedimentation...
 
          // remove particles sent left/right, coalesced or oud of domain and resize all n_part vectors
          if(opts.sedi || opts.adve || opts.coal)
            particles[dev_id].pimpl->hskpng_remove_n0();
 
          // updating particle->cell look-up table
          if (opts.adve || opts.sedi || n_rcyc)
            particles[dev_id].pimpl->hskpng_ijk();
        }
      }
      return res;
    }

    // diagnostic methods
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_sd_conc()
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].diag_sd_conc();
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].diag_dry_rng(r_mi, r_mx);
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].diag_wet_rng(r_mi, r_mx);
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_mom(const int &k)
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].diag_dry_mom(k);
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mom(const int &k)
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].diag_wet_mom(k);
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mass_dens(const real_t &a, const real_t &b)
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].diag_wet_mass_dens(a, b);
      }
    }

    // ...
//</listing>

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_chem(const enum chem_species_t &spec)
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].diag_chem(spec);
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_rw_ge_rc()
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].diag_rw_ge_rc();
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_RH_ge_Sc()
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].diag_RH_ge_Sc();
      }
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_all()
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].diag_all();
      }
    }

    template <typename real_t>
    real_t* particles_t<real_t, multi_CUDA>::outbuf()
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].pimpl->fill_outbuf();
        int n_cell_bfr;
        n_cell_bfr = dev_id * detail::get_dev_nx(glob_opts_init, 0) * detail::m1(glob_opts_init.ny) * detail::m1(glob_opts_init.nz);
        thrust::copy(
          particles[dev_id].pimpl->tmp_host_real_cell.begin(),
          particles[dev_id].pimpl->tmp_host_real_cell.end(),
          real_n_cell_tot.begin() + n_cell_bfr
        );
      }
      return &(*(real_n_cell_tot.begin()));
    }
  };
};
