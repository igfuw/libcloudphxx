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
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::xchng_domains(
    )
    {
#if defined(USE_MPI)
      // ranks of processes to the left/right, periodic boundary in x
      // TODO: same used in other places, unify
      const int lft_rank = mpi_rank > 0 ? mpi_rank - 1 : mpi_size - 1,
                rgt_rank = mpi_rank < mpi_size - 1 ? mpi_rank + 1 : 0;

      // exchange x0 x1 values
      for(int i=0; i<2; ++i)
      {
        // nonblocking send
        if( (i ? bcond.first : bcond.second) == detail::distmem_mpi)
        {
          printf("rank %d: i %d przed send\n", mpi_rank, i);
          MPI_Isend(
            i ? &opts_init.x0 : &opts_init.x1,
            1,                              // no of values
            detail::get_mpi_type<real_t>(), // type
            i ? lft_rank : rgt_rank,                       // dest comm
            i,                              // message tag
            MPI_COMM_WORLD,                 // communicator
            new MPI_Request()
          );
          printf("rank %d: i %d po send, x0 %lf x1 %lf\n", mpi_rank, i, opts_init.x0, opts_init.x1);
        }

          // blocking recv
        if( (i ? bcond.second : bcond.first) == detail::distmem_mpi)
        {
          printf("rank %d: i %d przed recv\n", mpi_rank, i);
          MPI_Recv(
            i ? &rgt_x0 : &lft_x1,
            1,                              // no of values
            detail::get_mpi_type<real_t>(), // type
            i ? rgt_rank : lft_rank,                       // src comm
            i,                              // message tag
            MPI_COMM_WORLD,                  // communicator
            MPI_STATUS_IGNORE
          );
          printf("rank %d: i %d po recv, lft_x1 %lf rgt_x0 %lf\n", mpi_rank, i, lft_x1, rgt_x0);
        }
//        printf("rank %d: i %d przed barrier\n", mpi_rank, i);
//        MPI_Barrier(MPI_COMM_WORLD);
//        printf("rank %d: i %d po barrier\n", mpi_rank, i);
      }
#endif
    }
  };
};
