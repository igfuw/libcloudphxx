// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */


// contains definitions of members of particles_t specialized for multiple GPUs
#include <omp.h>
#include "detail/remote_boundary.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    // time-stepping methods
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::step_sync(
      const opts_t<real_t> &opts,
      arrinfo_t<real_t> th,
      arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3,
      std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem
    )
    {
      #pragma omp parallel num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].step_sync(opts, th, rv, rhod, courant_1, courant_2, courant_3, ambient_chem);
      }
    }

    template <typename real_t>
    real_t particles_t<real_t, multi_CUDA>::step_async(
      const opts_t<real_t> &opts
    )
    {
      // cuda streams and events to control asynchronous copies
      // note: storing them in particles_multi_t caused errors
      // on program exit
      cudaStream_t streams[glob_opts_init.dev_count];
      cudaEvent_t events[glob_opts_init.dev_count];

      real_t res = 0.;
      #pragma omp parallel reduction(+:res) num_threads(glob_opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));

        // do step async on each device
        res = particles[dev_id].step_async(opts);

        // --- copy advected SDs to other devices ---
        if(opts.adve && glob_opts_init.dev_count>1)
        {
          namespace arg = thrust::placeholders;
          typedef unsigned long long n_t; // TODO: same typedef is in impl struct !! particles::impl::n_t ? 

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
          thrust_device::vector<real_t> &vt(particles[dev_id].pimpl->vt);
          thrust_device::vector<real_t> &out_real_bfr(particles[dev_id].pimpl->out_real_bfr);
          thrust_device::vector<real_t> &in_real_bfr(particles[dev_id].pimpl->in_real_bfr);
          thrust_device::vector<n_t> &n(particles[dev_id].pimpl->n);
          thrust_device::vector<n_t> &out_n_bfr(particles[dev_id].pimpl->out_n_bfr);
          thrust_device::vector<n_t> &in_n_bfr(particles[dev_id].pimpl->in_n_bfr);
          std::pair<detail::bcond_t, detai::bcond_t> &bcond(particles[dev_id].pimpl->bcond);
          // i and k must have not changed since impl->bcnd !!
          const thrust_device::vector<thrust_size_t> &lft_id(particles[dev_id].pimpl->i);
          const thrust_device::vector<thrust_size_t> &rgt_id(particles[dev_id].pimpl->k);

          // IDs of devices to the left/right, periodic_ext boundary in x
          const int lft_dev = dev_id > 0 ? dev_id - 1 : glob_opts_init.dev_count - 1,
                    rgt_dev = dev_id < glob_opts_init.dev_count-1 ? dev_id + 1 : 0;

          // init stream and event
          gpuErrchk(cudaStreamCreate(&streams[dev_id]));
          gpuErrchk(cudaEventCreateWithFlags(&events[dev_id], cudaEventDisableTiming ));

          // no of SDs copied
          int n_copied;
  
          // real_t vectors to be copied
          thrust_device::vector<real_t> * real_t_vctrs[] = {&rd3, &rw2, &kpa, &vt, &x, &z, &y};
          const int real_vctrs_count = glob_opts_init.ny == 0 ? 6 : 7;

          if(bcond.first == detail::distmem_cuda)
          {
            // prepare buffer with n_t to be copied left
            // TODO: serialize n_t and real_t with boost serialize
            thrust::copy(
              thrust::make_permutation_iterator(n.begin(), lft_id.begin()),
              thrust::make_permutation_iterator(n.begin(), lft_id.begin()) + lft_count,
              out_n_bfr.begin()
            );

            // start async copy of n buffer to the left
            gpuErrchk(cudaMemcpyPeerAsync(
              particles[lft_dev].pimpl->in_n_bfr.data().get(), lft_dev,  //dst
              out_n_bfr.data().get(), dev_id,                             //src 
              lft_count * sizeof(n_t),                                    //no of bytes
              streams[dev_id]                                             //best performance if stream belongs to src
            ));
            // record beginning of copying
            gpuErrchk(cudaEventRecord(events[dev_id], streams[dev_id]));
          }
          // barrier to make sure that all devices started copying
          #pragma omp barrier

          // adjust x of prtcls to be sent left to match new device's domain
          if(bcond.first == detail::distmem_cuda)
          {
            thrust::transform(
              thrust::make_permutation_iterator(x.begin(), lft_id.begin()),
              thrust::make_permutation_iterator(x.begin(), lft_id.begin()) + lft_count,
              thrust::make_permutation_iterator(x.begin(), lft_id.begin()), // in place
              detail::remote<real_t>(particles[dev_id].opts_init->x0, particles[lft_dev].opts_init->x1)
            );

            // prepare the real_t buffer for copy left
            for(int i = 0; i < real_vctrs_count; ++i)
              thrust::copy(
                thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), lft_id.begin()),
                thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), lft_id.begin()) + lft_count,
                out_real_bfr.begin() + i * lft_count
              );
          }

          if(bcond.second == detail::distmem_cuda)
          {
            // wait for the copy of n from right into current device to finish
            gpuErrchk(cudaEventSynchronize(events[rgt_dev]));

            // unpack the n buffer sent to this device from right
            n_copied = particles[rgt_dev].pimpl->lft_count;
            n_part_old = n_part;
            n_part += n_copied;
            n.resize(n_part);
            thrust::copy(in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);
          }

          // start async copy of real buffer to the left; same stream as n_bfr - will start only if previous copy finished
          if(bcond.first == detail::distmem_cuda)
          {
            gpuErrchk(cudaMemcpyPeerAsync(
              particles[lft_dev].pimpl->in_real_bfr.data().get(), lft_dev,  //dst
              out_real_bfr.data().get(), dev_id,                             //src 
              real_vctrs_count * lft_count * sizeof(real_t),                 //no of bytes
              streams[dev_id]                                                //best performance if stream belongs to src
            ));
            // record beginning of copying
            gpuErrchk(cudaEventRecord(events[dev_id], streams[dev_id]));
          }
          // barrier to make sure that all devices started copying
          #pragma omp barrier

          if(bcond.second == detail::distmem_cuda)
          {
            // prepare buffer with n_t to be copied right
            thrust::copy(
              thrust::make_permutation_iterator(n.begin(), rgt_id.begin()),
              thrust::make_permutation_iterator(n.begin(), rgt_id.begin()) + rgt_count,
              out_n_bfr.begin()
            );

            // adjust x of prtcls to be sent right to match new device's domain
            thrust::transform(
              thrust::make_permutation_iterator(x.begin(), rgt_id.begin()),
              thrust::make_permutation_iterator(x.begin(), rgt_id.begin()) + rgt_count,
              thrust::make_permutation_iterator(x.begin(), rgt_id.begin()), // in place
              detail::remote<real_t>(particles[dev_id].opts_init->x1, particles[rgt_dev].opts_init->x0)
            );

            // wait for the copy of real from right into current device to finish
            gpuErrchk(cudaEventSynchronize(events[rgt_dev]));

            // unpack the real buffer sent to this device from right
            for(int i = 0; i < real_vctrs_count; ++i)
            {
              real_t_vctrs[i]->resize(n_part);
              thrust::copy( in_real_bfr.begin() + i * n_copied, in_real_bfr.begin() + (i+1) * n_copied, real_t_vctrs[i]->begin() + n_part_old);
            }

            // start async copy of n buffer to the right
            gpuErrchk(cudaMemcpyPeerAsync(
              particles[rgt_dev].pimpl->in_n_bfr.data().get(), rgt_dev,  //dst
              out_n_bfr.data().get(), dev_id,                             //src 
              rgt_count * sizeof(n_t),                                    //no of bytes
              streams[dev_id]                                             //best performance if stream belongs to src
            ));
            // record beginning of copying
            gpuErrchk(cudaEventRecord(events[dev_id], streams[dev_id]));
            // barrier to make sure that all devices started copying
          }
          #pragma omp barrier

          // prepare the real_t buffer for copy to the right
          if(bcond.second == detail::distmem_cuda)
            for(int i = 0; i < real_vctrs_count; ++i)
              thrust::copy(
                thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), rgt_id.begin()),
                thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), rgt_id.begin()) + rgt_count,
                out_real_bfr.begin() + i * rgt_count
              );

          if(bcond.first == detail::distmem_cuda)
          {
            // wait for the copy of n from left into current device to finish
            gpuErrchk(cudaEventSynchronize(events[lft_dev]));
            // unpack the n buffer sent to this device from left
            n_copied = particles[lft_dev].pimpl->rgt_count;
            n_part_old = n_part;
            n_part += n_copied;
            n.resize(n_part);
            thrust::copy( in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);
          }

          if(bcond.second == detail::distmem_cuda)
          {
            // start async copy of real buffer to the right
            gpuErrchk(cudaMemcpyPeerAsync(
              particles[rgt_dev].pimpl->in_real_bfr.data().get(), rgt_dev,  //dst
              out_real_bfr.data().get(), dev_id,                             //src 
              real_vctrs_count * rgt_count * sizeof(real_t),                 //no of bytes
              streams[dev_id]                                                //best performance if stream belongs to src
            ));
            // record beginning of copying
            gpuErrchk(cudaEventRecord(events[dev_id], streams[dev_id]));
          }
          // barrier to make sure that all devices started copying
          #pragma omp barrier

          // flag SDs sent left/right for removal
          if(bcond.first == detail::distmem_cuda)
            thrust::copy(
              thrust::make_constant_iterator<n_t>(0),
              thrust::make_constant_iterator<n_t>(0) + lft_count,
              thrust::make_permutation_iterator(n.begin(), lft_id.begin())
            );
          if(bcond.second == detail::distmem_cuda)
            thrust::copy(
              thrust::make_constant_iterator<n_t>(0),
              thrust::make_constant_iterator<n_t>(0) + rgt_count,
              thrust::make_permutation_iterator(n.begin(), rgt_id.begin())
            );
          
          if(bcond.first == detail::distmem_cuda)
          {
            // wait for the copy of real from left into current device to finish
            gpuErrchk(cudaEventSynchronize(events[lft_dev]));

            // unpack the real buffer sent to this device from left
            for(int i = 0; i < real_vctrs_count; ++i)
            {
              real_t_vctrs[i]->resize(n_part);
              thrust::copy(in_real_bfr.begin() + i * n_copied, in_real_bfr.begin() + (i+1) * n_copied, real_t_vctrs[i]->begin() + n_part_old);
            }
          }
          // resize all vectors of size n_part
          particles[dev_id].pimpl->hskpng_resize_npart();

          // particles are not sorted now
          particles[dev_id].pimpl->sorted = false;          

          // clean streams and events
          #pragma omp barrier
          gpuErrchk(cudaStreamDestroy(streams[dev_id]));
          gpuErrchk(cudaEventDestroy(events[dev_id]));
        }

        // finalize async if more than one GPU on this node
        if(glob_opts_init.dev_count>1) 
          particles[dev_id].pimpl->post_copy(opts);
      }
      return res;
    }
  };
};
