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
    // fill courant number halos over mpi boundaries
    // TODO: similar to xchng_domains
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::xchng_courants(
    )
    {
#if defined(USE_MPI)
      if(halo_size == 0) return; // no halo to be exchanged

      // ranks of processes to the left/right, periodic boundary in x
      // TODO: same used in other places, unify
      const int lft_rank = mpi_rank > 0 ? mpi_rank - 1 : mpi_size - 1,
                rgt_rank = mpi_rank < mpi_size - 1 ? mpi_rank + 1 : 0;

      // indices of the locations from which courants should be copied
      const int cx_lft_internal_idx( 
          n_dims == 1 ? (halo_size + 1):                                      // 1D
          n_dims == 2 ? (halo_size + 1) * opts_init.nz:                       // 2D
                        (halo_size + 1) * opts_init.nz * opts_init.ny         // 3D
        );
      const int cx_rgt_internal_idx = n_cell; 

      const int cz_lft_internal_idx = halo_z; 
      const int cz_rgt_internal_idx( 
          n_dims == 2 ? opts_init.nx * (opts_init.nz + 1):                       // 2D
                        opts_init.nx * (opts_init.nz + 1) * opts_init.ny         // 3D
        );

      const int cy_lft_internal_idx = halo_y; 
      const int cy_rgt_internal_idx( 
                        opts_init.nx * opts_init.nz * (opts_init.ny + 1)         // 3D
        );

      // indices of the locations to which courants should be copied
      const int cx_lft_halo_idx = 0; 
      const int cx_rgt_halo_idx = courant_x.size() - halo_x; 
      const int cz_lft_halo_idx = 0; 
      const int cz_rgt_halo_idx = courant_z.size() - halo_z; 
      const int cy_lft_halo_idx = 0; 
      const int cy_rgt_halo_idx = courant_y.size() - halo_y; 

      enum { cx_tag = 12331, cy_tag, cz_tag };

      // exchange courant_X/Y/Z values
      for(int i=0; i<2; ++i) // left/right loop
      {
        // Courant x
        if(n_dims > 0)
        {
          // nonblocking send
          if( (i ? bcond.first : bcond.second) == detail::distmem_mpi)
          {
            MPI_Isend(
              i ? courant_x.data().get() + cx_lft_internal_idx : courant_x.data().get() +  cx_rgt_internal_idx,
              halo_x,                              // no of values
              detail::get_mpi_type<real_t>(),      // type
              i ? lft_rank : rgt_rank,             // dest comm
              (i+1) * cx_tag,                      // message tag
              detail::MPI_COMM_LIBCLOUD,           // communicator
              new MPI_Request()
            );
          }

          // blocking recv
          if( (i ? bcond.second : bcond.first) == detail::distmem_mpi)
          {
            MPI_Recv(
              i ? courant_x.data().get() + cx_rgt_halo_idx : courant_x.data().get() +  cx_lft_halo_idx,
              halo_x,                              // no of values
              detail::get_mpi_type<real_t>(), // type
              i ? rgt_rank : lft_rank,                       // src comm
              (i+1) * cx_tag,                      // message tag
              detail::MPI_COMM_LIBCLOUD,                  // communicator
              MPI_STATUS_IGNORE
            );
          }
        }
        // Courant z
        if(n_dims > 1)
        {
          // nonblocking send
          if( (i ? bcond.first : bcond.second) == detail::distmem_mpi)
          {
            MPI_Isend(
              i ? courant_z.data().get() + cz_lft_internal_idx : courant_z.data().get() +  cz_rgt_internal_idx,
              halo_z,                              // no of values
              detail::get_mpi_type<real_t>(),      // type
              i ? lft_rank : rgt_rank,             // dest comm
              (i+1) * cz_tag,                      // message tag
              detail::MPI_COMM_LIBCLOUD,           // communicator
              new MPI_Request()
            );
          }

          // blocking recv
          if( (i ? bcond.second : bcond.first) == detail::distmem_mpi)
          {
            MPI_Recv(
              i ? courant_z.data().get() + cz_rgt_halo_idx : courant_z.data().get() +  cz_lft_halo_idx,
              halo_z,                              // no of values
              detail::get_mpi_type<real_t>(), // type
              i ? rgt_rank : lft_rank,                       // src comm
              (i+1) * cz_tag,                      // message tag
              detail::MPI_COMM_LIBCLOUD,                  // communicator
              MPI_STATUS_IGNORE
            );
          }
        }
        // Courant y
        if(n_dims > 2)
        {
          // nonblocking send
          if( (i ? bcond.first : bcond.second) == detail::distmem_mpi)
          {
            MPI_Isend(
              i ? courant_y.data().get() + cy_lft_internal_idx : courant_y.data().get() +  cy_rgt_internal_idx,
              halo_y,                              // no of values
              detail::get_mpi_type<real_t>(),      // type
              i ? lft_rank : rgt_rank,             // dest comm
              (i+1) * cy_tag,                      // message tag
              detail::MPI_COMM_LIBCLOUD,           // communicator
              new MPI_Request()
            );
          }

          // blocking recv
          if( (i ? bcond.second : bcond.first) == detail::distmem_mpi)
          {
            MPI_Recv(
              i ? courant_y.data().get() + cy_rgt_halo_idx : courant_y.data().get() +  cy_lft_halo_idx,
              halo_y,                              // no of values
              detail::get_mpi_type<real_t>(), // type
              i ? rgt_rank : lft_rank,                       // src comm
              (i+1) * cy_tag,                      // message tag
              detail::MPI_COMM_LIBCLOUD,                  // communicator
              MPI_STATUS_IGNORE
            );
          }
        }
      }
#endif
    }
  };
};
