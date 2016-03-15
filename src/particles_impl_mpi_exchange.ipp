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
      enum {tag_n_t, tag_real_t};
      template <typename real_t>
      struct remote
      {
        real_t lcl, rmt;

        remote(real_t lcl, real_t rmt) : lcl(lcl), rmt(rmt) {}

        BOOST_GPU_ENABLED
        real_t operator()(real_t x)
        {
          real_t res = rmt + x - lcl;
          if(res == rmt) res = nextafter(res, real_t(0.)); // in single precision, we used to get x=x1
          return res;
        }
      };

      template<typename real_t>
      MPI_Datatype get_mpi_type()
      {
        throw std::runtime_error("Unsupported MPI datatype");
      }

      template<>
      MPI_datatype get_mpi_type<float>
      {
        return MPI_FLOAT;
      }

      template<>
      MPI_datatype get_mpi_type<double>
      {
        return MPI_DOUBLE;
      }

      template<>
      MPI_datatype get_mpi_type<unsigned long long>
      {
        return MPI_UNSIGNED_LONG_LONG;
      }
    };

    // --- copy advected SDs to other devices ---
    // TODO: many similarities to copy between GPUS in particles_impl_multi_gpu_step!
    template <typename real_t, typename backend_t>
    void particles_t<real_t, backend_t>::impl::mpi_exchange(
    )
    {
      namespace arg = thrust::placeholders;
      typedef unsigned long long n_t; // TODO: same typedef is in impl struct !! particles::impl::n_t ? 

      // i and k must have not changed since impl->bcnd !!
      const thrust_device::vector<thrust_size_t> &lft_id(i);
      const thrust_device::vector<thrust_size_t> &rgt_id(k);

      // ranks of processes to the left/right, periodic boundary in x
      const int lft_rank = mpi_rank > 0 ? mpi_rank - 1 : mpi_size - 1,
                rgt_rank = mpi_rank < mpi_size - 1 ? mpi_rank + 1 : 0;

      // prepare buffer with n_t to be copied left
      thrust::copy(
        thrust::make_permutation_iterator(n.begin(), lft_id.begin()),
        thrust::make_permutation_iterator(n.begin(), lft_id.begin()) + lft_count,
        out_n_bfr.begin()
      );

      // start async copy of n buffer to the left
      MPI_ISend(
        out_n_bfr.data().get(),       // raw pointer to the buffer
        lft_count,                    // no of values to send
        detail::get_mpi_type<n_t>,    // type
        lft_rank,                     // dest comm
        detail::tag_n_t,              // message tag
        MPI_COMM_WORLD                // communicator
      );

      // start async receiving of n buffer from right
      MPI_IRecv(
        in_n_bfr.data().get(),        // raw pointer to the buffer
        in_n_bfr.size(),              // max no of values to recv
        detail::get_mpi_type<n_t>,    // type
        rgt_rank,                     // src comm
        detail::tag_n_t,              // message tag
        MPI_COMM_WORLD,               // communicator
        MPI_STATUS_IGNORE             // status
      );


      // adjust x of prtcls to be sent left to match new device's domain
      thrust::transform(
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()),
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()) + lft_count,
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()), // in place
        detail::remote<real_t>(opts_init.x0, lft_dev_x1)
      );

      // prepare the real_t buffer for copy left
      thrust_device::vector<real_t> * real_t_vctrs[] = {&rd3, &rw2, &kpa, &vt, &x, &z, &y};
      const int real_vctrs_count = glob_opts_init.ny == 0 ? 6 : 7;
      for(int i = 0; i < real_vctrs_count; ++i)
        thrust::copy(
          thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), lft_id.begin()),
          thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), lft_id.begin()) + lft_count,
          out_real_bfr.begin() + i * lft_count
        );

      // start async copy of real buffer to the left
      MPI_ISend(
        out_real_bfr.data().get(),       // raw pointer to the buffer
        lft_count * real_vctrs_count,                    // no of values to send
        detail::get_mpi_type<real_t>,    // type
        lft_rank,                     // dest comm
        detail::tag_real_t,              // message tag
        MPI_COMM_WORLD                // communicator
      );

      // start async receiving of real buffer from right
      MPI_IRecv(
        in_real_bfr.data().get(),        // raw pointer to the buffer
        in_real_bfr.size(),              // max no of values to recv
        detail::get_mpi_type<real_t>,    // type
        rgt_rank,                     // src comm
        detail::tag_real_t,              // message tag
        MPI_COMM_WORLD,               // communicator
        MPI_STATUS_IGNORE             // status
      );

      // check if n buffer from right arrived
      MPI_Wait(...)

      // unpack the n buffer sent to this device from right
      int n_copied = particles[rgt_dev].pimpl->lft_count;
      n_part_old = n_part;
      n_part += n_copied;
      n.resize(n_part);
      thrust::copy(in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);

      // check if out_n_bfr sent left has been received

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
        detail::remote<real_t>(opts_init.x1, rgt_dev_x0)
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

      // resize all vectors of size n_part
      particles[dev_id].pimpl->hskpng_resize_npart();

      // particles are not sorted now
      particles[dev_id].pimpl->sorted = false;          

      // clean streams and events
      #pragma omp barrier
      gpuErrchk(cudaStreamDestroy(streams[dev_id]));
      gpuErrchk(cudaEventDestroy(events[dev_id]));
    }
  };
};
