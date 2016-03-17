// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include "detail/get_mpi_type.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      enum {tag_n_t = 654321, tag_real_t = 654322}; // hope other libs dont use these tags, TODO: using separate communicator would help?

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
    };

    // --- copy advected SDs to other devices ---
    // TODO: many similarities to copy between GPUS in particles_impl_multi_gpu_step!
    // TODO: use MPI's derived datatypes instead of packing/unpacking local buffers? wouldn't have to use separate buffers for n_t and real_t
    // TODO: use MPI's built-in [catresian] topology?
    // TODO: add MPI_CHECK over each send/recv/wait call
    // TODO: C++ MPI syntax not supported by all implementations? revert to C syntax?
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::mpi_exchange(
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

      // n_t/real_t send/receive requests
      MPI_Request req_recv_n_t, 
                  req_recv_real_t,
                  req_send_n_t, 
                  req_send_real_t;

      // mpi status handler
      MPI_Status status;

      // prepare buffer with n_t to be copied left
      thrust::copy(
        thrust::make_permutation_iterator(n.begin(), lft_id.begin()),
        thrust::make_permutation_iterator(n.begin(), lft_id.begin()) + lft_count,
        out_n_bfr.begin()
      );

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

      // start async receiving of n buffer from right
      MPI_Irecv(
        in_n_bfr.data().get(),        // raw pointer to the buffer
        in_n_bfr.size(),              // max no of values to recv
        detail::get_mpi_type<n_t>(),    // type
        rgt_rank,                     // src comm
        detail::tag_n_t,              // message tag
        MPI_COMM_WORLD,               // communicator
        &req_recv_n_t
      );


      // adjust x of prtcls to be sent left to match new device's domain
      thrust::transform(
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()),
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()) + lft_count,
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()), // in place
        detail::remote<real_t>(opts_init.x0, lft_x1)
      );

      // prepare the real_t buffer for copy left
      thrust_device::vector<real_t> * real_t_vctrs[] = {&rd3, &rw2, &kpa, &vt, &x, &z, &y};
      const int real_vctrs_count = opts_init.ny == 0 ? 6 : 7;
      for(int i = 0; i < real_vctrs_count; ++i)
        thrust::copy(
          thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), lft_id.begin()),
          thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), lft_id.begin()) + lft_count,
          out_real_bfr.begin() + i * lft_count
        );

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
      int n_copied;
      MPI_Get_count(&status, detail::get_mpi_type<n_t>(), &n_copied);
      n_part_old = n_part;
      n_part += n_copied;
      n.resize(n_part);
      thrust::copy(in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);

      // check if out_n_bfr sent left has been received
      MPI_Wait(&req_send_n_t, MPI_STATUS_IGNORE);

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
        detail::remote<real_t>(opts_init.x1, rgt_x0)
      );

      // wait for the copy of real from right into current device to finish
      MPI_Wait(&req_recv_real_t, MPI_STATUS_IGNORE);

      // unpack the real buffer sent to this device from right
      for(int i = 0; i < real_vctrs_count; ++i)
      {
        real_t_vctrs[i]->resize(n_part);
        thrust::copy( in_real_bfr.begin() + i * n_copied, in_real_bfr.begin() + (i+1) * n_copied, real_t_vctrs[i]->begin() + n_part_old);
      }

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

      // start async receiving of n buffer from left
      MPI_Irecv(
        in_n_bfr.data().get(),        // raw pointer to the buffer
        in_n_bfr.size(),              // max no of values to recv
        detail::get_mpi_type<n_t>(),    // type
        lft_rank,                     // src comm
        detail::tag_n_t,              // message tag
        MPI_COMM_WORLD,               // communicator
        &req_recv_n_t
      );

      // prepare the real_t buffer for copy to the right
      for(int i = 0; i < real_vctrs_count; ++i)
        thrust::copy(
          thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), rgt_id.begin()),
          thrust::make_permutation_iterator(real_t_vctrs[i]->begin(), rgt_id.begin()) + rgt_count,
          out_real_bfr.begin() + i * rgt_count
        );


      // check if n buffer from left arrived
      MPI_Wait(&req_recv_n_t, &status);

      // unpack the n buffer sent to this device from left
      MPI_Get_count(&status, detail::get_mpi_type<n_t>(), &n_copied);
      n_part_old = n_part;
      n_part += n_copied;
      n.resize(n_part);
      thrust::copy(in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);

      // start async copy of real buffer to the right
      MPI_Isend(
        out_real_bfr.data().get(),       // raw pointer to the buffer
        rgt_count * real_vctrs_count,                    // no of values to send
        detail::get_mpi_type<real_t>(),    // type
        rgt_rank,                     // dest comm
        detail::tag_real_t,              // message tag
        MPI_COMM_WORLD,                // communicator
        &req_send_real_t
      );

      // start async receiving of real buffer from left
      MPI_Irecv(
        in_real_bfr.data().get(),        // raw pointer to the buffer
        in_real_bfr.size(),              // max no of values to recv
        detail::get_mpi_type<real_t>(),    // type
        lft_rank,                     // src comm
        detail::tag_real_t,              // message tag
        MPI_COMM_WORLD,               // communicator
        &req_recv_real_t
      );

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
      MPI_Wait(&req_recv_real_t, MPI_STATUS_IGNORE);

      // unpack the real buffer sent to this device from left
      for(int i = 0; i < real_vctrs_count; ++i)
      {
        real_t_vctrs[i]->resize(n_part);
        thrust::copy(in_real_bfr.begin() + i * n_copied, in_real_bfr.begin() + (i+1) * n_copied, real_t_vctrs[i]->begin() + n_part_old);
      }

      // resize all vectors of size n_part
      hskpng_resize_npart();

      // particles are not sorted now
      sorted = false;          
    }
  };
};
