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
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].step_sync(opts, th, rv, courant_1, courant_2, courant_3, rhod);
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
          thrust_device::vector<real_t> &out_real_bfr(particles[dev_id].pimpl->out_real_bfr);
          thrust_device::vector<real_t> &in_real_bfr(particles[dev_id].pimpl->in_real_bfr);
          thrust_device::vector<n_t> &n(particles[dev_id].pimpl->n);
          thrust_device::vector<n_t> &out_n_bfr(particles[dev_id].pimpl->out_n_bfr);
          thrust_device::vector<n_t> &in_n_bfr(particles[dev_id].pimpl->in_n_bfr);
          // i and k must have not changed since impl->bcnd !!
          const thrust_device::vector<thrust_size_t> &lft_id(particles[dev_id].pimpl->i);
          const thrust_device::vector<thrust_size_t> &rgt_id(particles[dev_id].pimpl->k);

          // IDs of devices to the left/right, periodic boundary in x
          const int lft_dev = dev_id > 0 ? dev_id - 1 : glob_opts_init.dev_count - 1,
                    rgt_dev = dev_id < glob_opts_init.dev_count-1 ? dev_id + 1 : 0;

          // init stream and event
          gpuErrchk(cudaStreamCreate(&streams[dev_id]));
          gpuErrchk(cudaEventCreateWithFlags(&events[dev_id], cudaEventDisableTiming ));

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
          // barrier to make sure that all devices started copying
          #pragma omp barrier

          // adjust x of prtcls to be sent left to match new device's domain
          thrust::transform(
            thrust::make_permutation_iterator(x.begin(), lft_id.begin()),
            thrust::make_permutation_iterator(x.begin(), lft_id.begin()) + lft_count,
            thrust::make_permutation_iterator(x.begin(), lft_id.begin()), // in place
            arg::_1 + particles[lft_dev].opts_init->x1 - particles[dev_id].opts_init->x0   // operation
          );

          // prepare the real_t buffer for copy left
          thrust_device::vector<real_t> * real_t_vctrs[] = {&rd3, &rw2, &kpa, &x, &z, &y};
          const int real_vctrs_count = glob_opts_init.ny == 0 ? 5 : 6;
          for(int i = 0; i < real_vctrs_count; ++i)
            thrust::copy(
              thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), lft_id.begin()),
              thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), lft_id.begin()) + lft_count,
              out_real_bfr.begin() + i * lft_count
            );

          // wait for the copy of n from right into current device to finish
          gpuErrchk(cudaEventSynchronize(events[rgt_dev]));
          // unpack the n buffer sent to this device from right
          int n_copied = particles[rgt_dev].pimpl->lft_count;
          n_part_old = n_part;
          n_part += n_copied;
          n.resize(n_part);
          thrust::copy(in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);

          // start async copy of real buffer to the left; same stream as n_bfr - will start only if previous copy finished
          gpuErrchk(cudaMemcpyPeerAsync(
            particles[lft_dev].pimpl->in_real_bfr.data().get(), lft_dev,  //dst
            out_real_bfr.data().get(), dev_id,                             //src 
            real_vctrs_count * lft_count * sizeof(real_t),                 //no of bytes
            streams[dev_id]                                                //best performance if stream belongs to src
          ));
          // record beginning of copying
          gpuErrchk(cudaEventRecord(events[dev_id], streams[dev_id]));
          // barrier to make sure that all devices started copying
          #pragma omp barrier

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
            arg::_1 + particles[rgt_dev].opts_init->x0 - particles[dev_id].opts_init->x1   // operation
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
          #pragma omp barrier

          // prepare the real_t buffer for copy to the right
          for(int i = 0; i < real_vctrs_count; ++i)
            thrust::copy(
              thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), rgt_id.begin()),
              thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), rgt_id.begin()) + rgt_count,
              out_real_bfr.begin() + i * rgt_count
            );

          // wait for the copy of n from left into current device to finish
          gpuErrchk(cudaEventSynchronize(events[lft_dev]));
          // unpack the n buffer sent to this device from left
          n_copied = particles[lft_dev].pimpl->rgt_count;
          n_part_old = n_part;
          n_part += n_copied;
          n.resize(n_part);
          thrust::copy( in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);

          // start async copy of real buffer to the right
          gpuErrchk(cudaMemcpyPeerAsync(
            particles[rgt_dev].pimpl->in_real_bfr.data().get(), rgt_dev,  //dst
            out_real_bfr.data().get(), dev_id,                             //src 
            real_vctrs_count * rgt_count * sizeof(real_t),                 //no of bytes
            streams[dev_id]                                                //best performance if stream belongs to src
          ));
          // record beginning of copying
          gpuErrchk(cudaEventRecord(events[dev_id], streams[dev_id]));
          // barrier to make sure that all devices started copying
          #pragma omp barrier

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
          gpuErrchk(cudaEventSynchronize(events[lft_dev]));

          // unpack the real buffer sent to this device from left
          for(int i = 0; i < real_vctrs_count; ++i)
          {
            real_t_vctrs[i]->resize(n_part);
            thrust::copy(in_real_bfr.begin() + i * n_copied, in_real_bfr.begin() + (i+1) * n_copied, real_t_vctrs[i]->begin() + n_part_old);
          }

          // particles are not sorted now
          particles[dev_id].pimpl->sorted = false;          

          // clean streams and events
          #pragma omp barrier
          gpuErrchk(cudaStreamDestroy(streams[dev_id]));
          gpuErrchk(cudaEventDestroy(events[dev_id]));
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
  };
};
