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
    // --- copy advected SDs to other devices ---
    // TODO: many similarities to copy between GPUS in particles_impl_multi_gpu_step!
    // TODO: use MPI's derived datatypes instead of packing/unpacking local buffers? wouldn't have to use separate buffers for n_t and real_t
    // TODO: use MPI's built-in [catresian] topology?
    // TODO: add MPI_CHECK over each send/recv/wait call
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::mpi_exchange(
    )
    {
#if defined(USE_MPI)
      namespace arg = thrust::placeholders;

      // ranks of processes to the left/right, periodic boundary in x
      const int lft_rank = mpi_rank > 0 ? mpi_rank - 1 : mpi_size - 1,
                rgt_rank = mpi_rank < mpi_size - 1 ? mpi_rank + 1 : 0;

      // n_t/real_t send/receive requests
      MPI_Request req_recv_n_t, 
                  req_recv_real_t,
                  req_send_n_t, 
                  req_send_real_t;

      // mpi status handler
      MPI_Status status;

      if(bcond.first == detail::distmem_mpi)
      {
        // prepare buffer with n_t to be copied left
        pack_n_lft();

        // start async copy of n buffer to the left
        MPI_Isend(
          out_n_bfr.data().get(),       // raw pointer to the buffer
          lft_count,                    // no of values to send
          detail::get_mpi_type<n_t>(),    // type
          lft_rank,                     // dest comm
          detail::tag_n_t,              // message tag
          MPI_COMM_WORLD,
          &req_send_n_t
        );
      }
      // start async receiving of n buffer from right
      if(bcond.second == detail::distmem_mpi)
      {
        MPI_Irecv(
          in_n_bfr.data().get(),        // raw pointer to the buffer
          in_n_bfr.size(),              // max no of values to recv
          detail::get_mpi_type<n_t>(),    // type
          rgt_rank,                     // src comm
          detail::tag_n_t,              // message tag
          MPI_COMM_WORLD,               // communicator
          &req_recv_n_t
        );
      }

      if(bcond.first == detail::distmem_mpi)
      {
        // adjust x of prtcls to be sent left to match new device's domain
        bcnd_remote_lft(opts_init.x0, lft_x1);

        // prepare the real_t buffer for copy left
        pack_real_lft();

        // start async copy of real buffer to the left
        MPI_Isend(
          out_real_bfr.data().get(),       // raw pointer to the buffer
          lft_count * real_vctrs_count,                    // no of values to send
          detail::get_mpi_type<real_t>(),    // type
          lft_rank,                     // dest comm
          detail::tag_real_t,              // message tag
          MPI_COMM_WORLD,                // communicator
          &req_send_real_t
        );
      }

      if(bcond.second == detail::distmem_mpi)
      {
        // start async receiving of real buffer from right
        MPI_Irecv(
          in_real_bfr.data().get(),        // raw pointer to the buffer
          in_real_bfr.size(),              // max no of values to recv
          detail::get_mpi_type<real_t>(),    // type
          rgt_rank,                     // src comm
          detail::tag_real_t,              // message tag
          MPI_COMM_WORLD,               // communicator
          &req_recv_real_t
        );

        // check if n buffer from right arrived
        MPI_Wait(&req_recv_n_t, &status);

        // unpack the n buffer sent to this device from right
        MPI_Get_count(&status, detail::get_mpi_type<n_t>(), &n_copied);
        unpack_n(n_copied);
      }

      // check if out_n_bfr sent left has been received
      if(bcond.first == detail::distmem_mpi)
      {
        MPI_Wait(&req_send_n_t, MPI_STATUS_IGNORE);
      }

      if(bcond.second == detail::distmem_mpi)
      {
        // prepare buffer with n_t to be copied right
        pack_n_rgt();

        // adjust x of prtcls to be sent right to match new device's domain
        bcnd_remote_rgt(opts_init.x1, rgt_x0);

        // wait for the copy of real from right into current device to finish
        MPI_Wait(&req_recv_real_t, MPI_STATUS_IGNORE);

        // unpack the real buffer sent to this device from right
        unpack_real(n_copied);

        // start async copy of n buffer to the right
        MPI_Isend(
          out_n_bfr.data().get(),       // raw pointer to the buffer
          rgt_count,                    // no of values to send
          detail::get_mpi_type<n_t>(),    // type
          rgt_rank,                     // dest comm
          detail::tag_n_t,              // message tag
          MPI_COMM_WORLD,                // communicator
          &req_send_n_t
        );
      }

      // start async receiving of n buffer from left
      if(bcond.first == detail::distmem_mpi)
      {
        MPI_Irecv(
          in_n_bfr.data().get(),        // raw pointer to the buffer
          in_n_bfr.size(),              // max no of values to recv
          detail::get_mpi_type<n_t>(),    // type
          lft_rank,                     // src comm
          detail::tag_n_t,              // message tag
          MPI_COMM_WORLD,               // communicator
          &req_recv_n_t
        );
      }

      // prepare the real_t buffer for copy to the right
      if(bcond.second == detail::distmem_mpi)
        pack_real_rgt();

      // check if n buffer from left arrived
      if(bcond.first == detail::distmem_mpi)
      {
        MPI_Wait(&req_recv_n_t, &status);

        // unpack the n buffer sent to this device from left
        MPI_Get_count(&status, detail::get_mpi_type<n_t>(), &n_copied);
        unpack_n(n_copied);
      }

      // start async copy of real buffer to the right
      if(bcond.second == detail::distmem_mpi)
      {
        MPI_Isend(
          out_real_bfr.data().get(),       // raw pointer to the buffer
          rgt_count * real_vctrs_count,                    // no of values to send
          detail::get_mpi_type<real_t>(),    // type
          rgt_rank,                     // dest comm
          detail::tag_real_t,              // message tag
          MPI_COMM_WORLD,                // communicator
          &req_send_real_t
        );
      }

      // start async receiving of real buffer from left
      if(bcond.first == detail::distmem_mpi)
      {
        MPI_Irecv(
          in_real_bfr.data().get(),        // raw pointer to the buffer
          in_real_bfr.size(),              // max no of values to recv
          detail::get_mpi_type<real_t>(),    // type
          lft_rank,                     // src comm
          detail::tag_real_t,              // message tag
          MPI_COMM_WORLD,               // communicator
          &req_recv_real_t
        );
      }

      // flag SDs sent left/right for removal
      if(bcond.first == detail::distmem_mpi)
        flag_lft();

      if(bcond.second == detail::distmem_mpi)
        flag_rgt();
      

      // wait for the copy of real from left into current device to finish
      if(bcond.first == detail::distmem_mpi)
      {
        MPI_Wait(&req_recv_real_t, MPI_STATUS_IGNORE);

        // unpack the real buffer sent to this device from left
        unpack_real(n_copied);
      }

      // resize all vectors of size n_part
      hskpng_resize_npart();

      // particles are not sorted now
      sorted = false;          

      // wait for all  processes to finish before continuing
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
  };
};
