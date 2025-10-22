// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <chrono>
#include <thread>

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
//      MPI_CHECK(MPI_Barrier(detail::MPI_COMM_LIBCLOUD));

      if(!distmem_mpi()) return;
      namespace arg = thrust::placeholders;

      // ranks of processes to the left/right, periodic boundary in x
      const int lft_rank = mpi_rank > 0 ? mpi_rank - 1 : mpi_size - 1,
                rgt_rank = mpi_rank < mpi_size - 1 ? mpi_rank + 1 : 0;

      // n_t/real_t send/receive requests
      MPI_Request req_recv_n_t, 
                  req_recv_real_t,
                  req_send_n_t, 
                  req_send_real_t;

      // no of SDs copied
      int n_copied;

      // mpi status handler
      MPI_Status status;

//      MPI_CHECK(MPI_Barrier(detail::MPI_COMM_LIBCLOUD));

      // debugging information
//#if !defined(NDEBUG)
//      std::vector<real_t> dbg_sum;
//      dbg_sum.resize(distmem_real_vctrs.size());
//
//      MPI_Request req_dbg_recv_n_t, 
//                  req_dbg_recv_real_t,
//                  req_dbg_send_n_t, 
//                  req_dbg_send_real_t;
//#endif


      if(bcond.first == detail::distmem_mpi)
      {
        // prepare buffer with n_t to be copied left
        pack_n_lft();

        // without synchronize, we sometimes get errors in the MPI copy. E.g. rd3 gets too big (assert in cond_common finds this; tested on dycoms short test on Prometheus)
        // Is this because thrust::copy in pack_n_lft is not synchronous? 
        // At least pre Thrust 1.9.4. 1.9.4 made all calls (that are not explicitly asynchronous) synchronous, see https://github.com/thrust/thrust/blob/master/doc/changelog.md#thrust-194-cuda-101
        // TODO: Test with Thurst >= 1.9.4 if this sync is still needed
#if defined(__NVCC__)
        gpuErrchk(cudaDeviceSynchronize());
#endif

//        std::cerr << "mpi exchange: sending n left, lft_count = " << lft_count << " sum of out_n_bfr = " << thrust::reduce(out_n_bfr.begin(), out_n_bfr.begin() + lft_count) << std::endl;

        // start async copy of n buffer to the left
        MPI_CHECK(MPI_Isend(
          out_n_bfr.data().get(),       // raw pointer to the buffer
          lft_count * distmem_n_vctrs.size(),                    // no of values to send
          detail::get_mpi_type<n_t>(),    // type
          lft_rank,                     // dest comm
          detail::tag_n_lft,              // message tag
          detail::MPI_COMM_LIBCLOUD,
          &req_send_n_t
        ));
      }
//      MPI_CHECK(MPI_Barrier(detail::MPI_COMM_LIBCLOUD));
      // start async receiving of n buffer from right
      if(bcond.second == detail::distmem_mpi)
      {
        MPI_CHECK(MPI_Irecv(
          in_n_bfr.data().get(),        // raw pointer to the buffer
          in_n_bfr.size(),              // max no of values to recv
          detail::get_mpi_type<n_t>(),    // type
          rgt_rank,                     // src comm
          detail::tag_n_lft,              // message tag
          detail::MPI_COMM_LIBCLOUD,               // communicator
          &req_recv_n_t
        ));
      }

      if(bcond.first == detail::distmem_mpi)
      {
        // adjust x of prtcls to be sent left to match new device's domain
        bcnd_remote_lft(opts_init.x0, lft_x1);

        // prepare the real_t buffer for copy left
        pack_real_lft();

//#if !defined(NDEBUG)
//      // calculate sum of each vector for debugging
//      auto it = distmem_real_vctrs.begin();
//      while (it != distmem_real_vctrs.end())
//      {
//        auto dist = std::distance(distmem_real_vctrs.begin(), it);
//        dbg_sum.at(dist) = 
//          thrust::reduce(
//            out_real_bfr.begin() + dist * lft_count,
//            out_real_bfr.begin() + (dist+1) * lft_count
//          );
//
//        MPI_CHECK(MPI_Isend(
//          dbg_sum.data(),       // raw pointer to the buffer
//          distmem_real_vctrs.size(),                    // no of values to send
//          detail::get_mpi_type<real_t>(),    // type
//          lft_rank,                     // dest comm
//          detail::tag_n_lft,              // message tag
//          detail::MPI_COMM_LIBCLOUD,
//          &req_dbg_send_real_t
//        ));
//      }
//#endif

//      check comment for the previous device sync
#if defined(__NVCC__)
        gpuErrchk(cudaDeviceSynchronize());
#endif

//        std::cerr << "mpi exchange: sending real left, lft_count = " << lft_count << " sum of out_real_bfr = " << thrust::reduce(out_real_bfr.begin(), out_real_bfr.begin() + lft_count * distmem_real_vctrs.size()) << std::endl;

        // start async copy of real buffer to the left
        MPI_CHECK(MPI_Isend(
          out_real_bfr.data().get(),       // raw pointer to the buffer
          lft_count * distmem_real_vctrs.size(),                    // no of values to send
          detail::get_mpi_type<real_t>(),    // type
          lft_rank,                     // dest comm
          detail::tag_real_lft,              // message tag
          detail::MPI_COMM_LIBCLOUD,                // communicator
          &req_send_real_t
        ));
      }
//      MPI_CHECK(MPI_Barrier(detail::MPI_COMM_LIBCLOUD));

      if(bcond.second == detail::distmem_mpi)
      {
        // start async receiving of real buffer from right
        MPI_CHECK(MPI_Irecv(
          in_real_bfr.data().get(),        // raw pointer to the buffer
          in_real_bfr.size(),              // max no of values to recv
          detail::get_mpi_type<real_t>(),    // type
          rgt_rank,                     // src comm
          detail::tag_real_lft,              // message tag
          detail::MPI_COMM_LIBCLOUD,               // communicator
          &req_recv_real_t
        ));

        // check if n buffer from right arrived
        MPI_CHECK(MPI_Wait(&req_recv_n_t, &status));

        // unpack the n buffer sent to this device from right
        MPI_CHECK(MPI_Get_count(&status, detail::get_mpi_type<n_t>(), &n_copied));

//        std::cerr << "mpi exchange: receiving n from right, n_copied = " << n_copied << " sum of in_n_bfr = " << thrust::reduce(in_n_bfr.begin(), in_n_bfr.begin() + n_copied) << std::endl;
        unpack_n(n_copied);
      }

      if(bcond.second == detail::distmem_mpi)
      {
        // check if out_n_bfr sent left has been received
        if(bcond.first == detail::distmem_mpi)
          MPI_CHECK(MPI_Wait(&req_send_n_t, MPI_STATUS_IGNORE));

//        std::this_thread::sleep_for(std::chrono::seconds(1));

        // prepare buffer with n_t to be copied right
        pack_n_rgt();
        // no cudaDeviceSynchronize after this pack, because there's plenty of calls before ISend...

        // adjust x of prtcls to be sent right to match new device's domain
        bcnd_remote_rgt(opts_init.x1, rgt_x0);

        // wait for the copy of real from right into current device to finish
        MPI_CHECK(MPI_Wait(&req_recv_real_t, MPI_STATUS_IGNORE));
//        MPI_CHECK(MPI_Wait(&req_recv_real_t, &status));
  //      MPI_CHECK(MPI_Get_count(&status, detail::get_mpi_type<real_t>(), &n_copied));


//        std::cerr << "mpi exchange: receiving real from right, n_copied = " << n_copied << " sum of in_real_bfr = " << thrust::reduce(in_real_bfr.begin(), in_real_bfr.begin() + n_copied * distmem_real_vctrs.size()) << std::endl;

        // unpack the real buffer sent to this device from right
        unpack_real(n_copied);

//        std::cerr << "mpi exchange: sending n rgt, rgt_count = " << rgt_count <<  " sum of out_n_bfr = " << thrust::reduce(out_n_bfr.begin(), out_n_bfr.begin() + rgt_count) << std::endl;

/*
        std::cerr << "mpi exchange: sending n rgt, rgt_count = " << rgt_count << std::endl 
                  << " out_n_bfr.data.get = " << out_n_bfr.data().get() << std::endl
                  << " val at out_n_bfr.data.get = " << *out_n_bfr.data() << std::endl
                  << " out_n_bfr.data.get + rgt_count = " << out_n_bfr.data().get() + rgt_count << std::endl
                  << " val at (out_n_bfr.data.get + rgt_count) = " << *(out_n_bfr.data() + rgt_count) << std::endl
                  << " rgt_rank = " << rgt_rank << std::endl;
        //std::this_thread::sleep_for(std::chrono::milliseconds(10000));
        std::cerr << "mpi exchange after sleep: sending n rgt, rgt_count = " << rgt_count << std::endl 
                  << " out_n_bfr.data.get = " << out_n_bfr.data().get() << std::endl
                  << " val at out_n_bfr.data.get = " << *out_n_bfr.data() << std::endl
                  << " out_n_bfr.data.get + rgt_count = " << out_n_bfr.data().get() + rgt_count << std::endl
                  << " val at (out_n_bfr.data.get + rgt_count) = " << *(out_n_bfr.data() + rgt_count) << std::endl
                  << " rgt_rank = " << rgt_rank << std::endl;
*/
        // start async copy of n buffer to the right
        MPI_CHECK(MPI_Isend(
          out_n_bfr.data().get(),       // raw pointer to the buffer
          rgt_count * distmem_n_vctrs.size(),                    // no of values to send
          detail::get_mpi_type<n_t>(),    // type
          rgt_rank,                     // dest comm
          detail::tag_n_rgt,              // message tag
          detail::MPI_COMM_LIBCLOUD,                // communicator
          &req_send_n_t
        ));
      }
//      MPI_CHECK(MPI_Barrier(detail::MPI_COMM_LIBCLOUD));

      // start async receiving of n buffer from left
      if(bcond.first == detail::distmem_mpi)
      {
        MPI_CHECK(MPI_Irecv(
          in_n_bfr.data().get(),        // raw pointer to the buffer
          in_n_bfr.size(),              // max no of values to recv
          detail::get_mpi_type<n_t>(),    // type
          lft_rank,                     // src comm
          detail::tag_n_rgt,              // message tag
          detail::MPI_COMM_LIBCLOUD,               // communicator
          &req_recv_n_t
        ));
      }

      // prepare the real_t buffer for copy to the right
      if(bcond.second == detail::distmem_mpi)
      {
        // check if real_t buffer sent left has been received
        if(bcond.first == detail::distmem_mpi)
        {
          MPI_CHECK(MPI_Wait(&req_send_real_t, MPI_STATUS_IGNORE));
        }
        pack_real_rgt();
        // no cudaDeviceSynchronize after this pack, because there's plenty of calls before ISend...
      }

      // check if n buffer from left arrived
      if(bcond.first == detail::distmem_mpi)
      {
        MPI_CHECK(MPI_Wait(&req_recv_n_t, &status));

        // unpack the n buffer sent to this device from left
        MPI_CHECK(MPI_Get_count(&status, detail::get_mpi_type<n_t>(), &n_copied));

//        std::cerr << "mpi exchange: receiving n from left, n_copied = " << n_copied << " sum of in_n_bfr = " << thrust::reduce(in_n_bfr.begin(), in_n_bfr.begin() + n_copied) << std::endl;
        unpack_n(n_copied);
      }

      // start async copy of real buffer to the right
      if(bcond.second == detail::distmem_mpi)
      {

//        std::cerr << "mpi exchange: sending real rgt, rgt_count = " << rgt_count << " sum of out_real_bfr = " << thrust::reduce(out_real_bfr.begin(), out_real_bfr.begin() + rgt_count * distmem_real_vctrs.size()) << std::endl;

        MPI_CHECK(MPI_Isend(
          out_real_bfr.data().get(),       // raw pointer to the buffer
          rgt_count * distmem_real_vctrs.size(),                    // no of values to send
          detail::get_mpi_type<real_t>(),    // type
          rgt_rank,                     // dest comm
          detail::tag_real_rgt,              // message tag
          detail::MPI_COMM_LIBCLOUD,                // communicator
          &req_send_real_t
        ));
      }
//      MPI_CHECK(MPI_Barrier(detail::MPI_COMM_LIBCLOUD));

      // start async receiving of real buffer from left
      if(bcond.first == detail::distmem_mpi)
      {
        MPI_CHECK(MPI_Irecv(
          in_real_bfr.data().get(),        // raw pointer to the buffer
          in_real_bfr.size(),              // max no of values to recv
          detail::get_mpi_type<real_t>(),    // type
          lft_rank,                     // src comm
          detail::tag_real_rgt,              // message tag
          detail::MPI_COMM_LIBCLOUD,               // communicator
          &req_recv_real_t
        ));
      }

      // flag SDs sent left/right for removal
      if(bcond.first == detail::distmem_mpi)
        flag_lft();

      if(bcond.second == detail::distmem_mpi)
        flag_rgt();

      // wait for the copy of real from left into current device to finish
      if(bcond.first == detail::distmem_mpi)
      {
        MPI_CHECK(MPI_Wait(&req_recv_real_t, MPI_STATUS_IGNORE));

        // unpack the real buffer sent to this device from left

//        std::cerr << "mpi exchange: receiving real from left, n_copied = " << n_copied << " sum of in_real_bfr = " << thrust::reduce(in_real_bfr.begin(), in_real_bfr.begin() + n_copied * distmem_real_vctrs.size()) << std::endl;
        unpack_real(n_copied);
      }

      // resize all vectors of size n_part
      hskpng_resize_npart();

      // particles are not sorted now
      sorted = false;          

      // wait for all sends to finish to avoid external overwriting of the send buffer (e.g. multi_CUDA intra-node communications)
      MPI_CHECK(MPI_Wait(&req_send_n_t, MPI_STATUS_IGNORE));
      MPI_CHECK(MPI_Wait(&req_send_real_t, MPI_STATUS_IGNORE));
      //MPI_CHECK(MPI_Barrier(detail::MPI_COMM_LIBCLOUD));
      
#endif
    }
  };
};
