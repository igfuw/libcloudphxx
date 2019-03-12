// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct nextafter_fctr
      {
        real_t goal;
        nextafter_fctr(real_t goal): goal(goal) {}
        BOOST_GPU_ENABLED
        real_t operator()(real_t x)
        {
          return nextafter(x, goal);
        }
      };
    };

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::impl::step_async_and_copy(
      const opts_t<real_t> &opts,
      const int dev_id,
      std::vector<cudaStream_t> &streams,
      std::vector<cudaEvent_t> &events,
      detail::barrier_t &barrier
    )
    {
      gpuErrchk(cudaSetDevice(dev_id));

      // do step async on each device
      particles[dev_id]->step_async(opts);

      // --- copy advected SDs to other devices on the same node ---
      if(opts.adve && glob_opts_init.dev_count>1)
      {
        namespace arg = thrust::placeholders;
        typedef unsigned long long n_t; // TODO: same typedef is in impl struct !! particles::impl::n_t ? 

        // helper aliases
        const thrust_size_t &lft_count(particles[dev_id]->pimpl->lft_count);
        const thrust_size_t &rgt_count(particles[dev_id]->pimpl->rgt_count);
        auto &n_part(particles[dev_id]->pimpl->n_part);
        auto &n_part_old(particles[dev_id]->pimpl->n_part_old);
 	const int &distmem_real_vctrs_count(particles[dev_id]->pimpl->distmem_real_vctrs.size());
        thrust_device::vector<real_t> &sstp_tmp_th(particles[dev_id]->pimpl->sstp_tmp_th);
        thrust_device::vector<real_t> &sstp_tmp_rh(particles[dev_id]->pimpl->sstp_tmp_rh);
        thrust_device::vector<real_t> &sstp_tmp_rv(particles[dev_id]->pimpl->sstp_tmp_rv);
        thrust_device::vector<real_t> &sstp_tmp_p(particles[dev_id]->pimpl->sstp_tmp_p);
        thrust_device::vector<real_t> &out_real_bfr(particles[dev_id]->pimpl->out_real_bfr);
        thrust_device::vector<real_t> &in_real_bfr(particles[dev_id]->pimpl->in_real_bfr);
        thrust_device::vector<real_t> &x(particles[dev_id]->pimpl->x);
        thrust_device::vector<n_t> &out_n_bfr(particles[dev_id]->pimpl->out_n_bfr);
        thrust_device::vector<n_t> &in_n_bfr(particles[dev_id]->pimpl->in_n_bfr);
	std::pair<detail::bcond_t, detail::bcond_t> &bcond(particles[dev_id]->pimpl->bcond);

        // IDs of devices to the left/right, periodic_ext boundary in x
        const int lft_dev = dev_id > 0 ? dev_id - 1 : glob_opts_init.dev_count - 1,
                  rgt_dev = dev_id < glob_opts_init.dev_count-1 ? dev_id + 1 : 0;

        // init stream and event
        gpuErrchk(cudaStreamCreate(&streams[dev_id]));
        gpuErrchk(cudaEventCreateWithFlags(&events[dev_id], cudaEventDisableTiming ));

        // prepare buffer with n_t to be copied left
        if(bcond.first == detail::distmem_cuda)
        {
          // prepare buffer with n_t to be copied left
          // TODO: serialize n_t and real_t with boost serialize
          particles[dev_id]->pimpl->pack_n_lft();

          // start async copy of n buffer to the left
          gpuErrchk(cudaMemcpyPeerAsync(
            particles[lft_dev]->pimpl->in_n_bfr.data().get(), lft_dev,  //dst
            out_n_bfr.data().get(), dev_id,                             //src 
            lft_count * sizeof(n_t),                                    //no of bytes
            streams[dev_id]                                             //best performance if stream belongs to src
          ));
          // record beginning of copying
          gpuErrchk(cudaEventRecord(events[dev_id], streams[dev_id]));
        }

        // barrier to make sure that all devices started copying
        barrier.wait();

        // adjust x of prtcls to be sent left to match new device's domain
        if(bcond.first == detail::distmem_cuda)
        {
          // adjust x of prtcls to be sent left to match new device's domain
          particles[dev_id]->pimpl->bcnd_remote_lft(particles[dev_id]->opts_init->x0, particles[lft_dev]->opts_init->x1);

          // prepare the real_t buffer for copy left
          particles[dev_id]->pimpl->pack_real_lft();
        }

        if(bcond.second == detail::distmem_cuda)
        {
          // wait for the copy of n from right into current device to finish
          gpuErrchk(cudaEventSynchronize(events[rgt_dev]));
          // unpack the n buffer sent to this device from right
          particles[dev_id]->pimpl->unpack_n(particles[rgt_dev]->pimpl->lft_count); // also sets n_part_old and n_part
        }

        // start async copy of real buffer to the left; same stream as n_bfr - will start only if previous copy finished
        if(bcond.first == detail::distmem_cuda)
        {
          gpuErrchk(cudaMemcpyPeerAsync(
            particles[lft_dev]->pimpl->in_real_bfr.data().get(), lft_dev,  //dst
            out_real_bfr.data().get(), dev_id,                             //src 
            distmem_real_vctrs_count * lft_count * sizeof(real_t),                 //no of bytes
            streams[dev_id]                                                //best performance if stream belongs to src
          ));
          // record beginning of copying
          gpuErrchk(cudaEventRecord(events[dev_id], streams[dev_id]));
        }
        // barrier to make sure that all devices started copying
        barrier.wait();

        if(bcond.second == detail::distmem_cuda)
        {
          // prepare buffer with n_t to be copied right
          particles[dev_id]->pimpl->pack_n_rgt();

          // adjust x of prtcls to be sent right to match new device's domain
          particles[dev_id]->pimpl->bcnd_remote_rgt(particles[dev_id]->opts_init->x1, particles[rgt_dev]->opts_init->x0);

          // wait for the copy of real from right into current device to finish
          gpuErrchk(cudaEventSynchronize(events[rgt_dev]));

          // unpack the real buffer sent to this device from right
          particles[dev_id]->pimpl->unpack_real(particles[rgt_dev]->pimpl->lft_count);

          // sanitize x==x1 that could happen due to errors in copying?
          thrust::transform_if(x.begin() + n_part_old, x.begin() + n_part, x.begin() + n_part_old, detail::nextafter_fctr<real_t>(0.), arg::_1 == particles[dev_id]->opts_init->x1);

          // start async copy of n buffer to the right
          gpuErrchk(cudaMemcpyPeerAsync(
            particles[rgt_dev]->pimpl->in_n_bfr.data().get(), rgt_dev,  //dst
            out_n_bfr.data().get(), dev_id,                             //src 
            rgt_count * sizeof(n_t),                                    //no of bytes
            streams[dev_id]                                             //best performance if stream belongs to src
          ));
          // record beginning of copying
          gpuErrchk(cudaEventRecord(events[dev_id], streams[dev_id]));
          // barrier to make sure that all devices started copying
        }
        barrier.wait();

        // prepare the real_t buffer for copy to the right
        if(bcond.second == detail::distmem_cuda)
           particles[dev_id]->pimpl->pack_real_rgt();

        if(bcond.first == detail::distmem_cuda)
        {
          // wait for the copy of n from left into current device to finish
          gpuErrchk(cudaEventSynchronize(events[lft_dev]));
          // unpack the n buffer sent to this device from left
          particles[dev_id]->pimpl->unpack_n(particles[lft_dev]->pimpl->rgt_count); // also sets n_part etc..
        }

        if(bcond.second == detail::distmem_cuda)
        {
          // start async copy of real buffer to the right
          gpuErrchk(cudaMemcpyPeerAsync(
            particles[rgt_dev]->pimpl->in_real_bfr.data().get(), rgt_dev,  //dst
            out_real_bfr.data().get(), dev_id,                             //src 
            distmem_real_vctrs_count * rgt_count * sizeof(real_t),                 //no of bytes
            streams[dev_id]                                                //best performance if stream belongs to src
          ));
          // record beginning of copying
          gpuErrchk(cudaEventRecord(events[dev_id], streams[dev_id]));
        }

        // barrier to make sure that all devices started copying
        barrier.wait();
        if(bcond.first == detail::distmem_cuda)
          particles[dev_id]->pimpl->flag_lft(); 
        if(bcond.second == detail::distmem_cuda)
          particles[dev_id]->pimpl->flag_rgt(); 
        
        if(bcond.first == detail::distmem_cuda)
        {
          // wait for the copy of real from left into current device to finish
          gpuErrchk(cudaEventSynchronize(events[lft_dev]));

          // unpack the real buffer sent to this device from left
          particles[dev_id]->pimpl->unpack_real(particles[lft_dev]->pimpl->rgt_count);
        }

        // resize all vectors of size n_part
        particles[dev_id]->pimpl->hskpng_resize_npart();

        // particles are not sorted now
        particles[dev_id]->pimpl->sorted = false;          

        // clean streams and events
        barrier.wait();
        gpuErrchk(cudaStreamDestroy(streams[dev_id]));
        gpuErrchk(cudaEventDestroy(events[dev_id]));
      }
      // finalize async
      particles[dev_id]->pimpl->post_copy(opts);
    }
  };
};
