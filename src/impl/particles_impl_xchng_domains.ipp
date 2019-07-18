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
          MPI_Isend(
            i ? &opts_init.x0 : &opts_init.x1,
            1,                              // no of values
            detail::get_mpi_type<real_t>(), // type
            i ? lft_rank : rgt_rank,                       // dest comm
            i,                              // message tag
            detail::MPI_COMM_LIBCLOUD,                 // communicator
            new MPI_Request()
          );
        }

          // blocking recv
        if( (i ? bcond.second : bcond.first) == detail::distmem_mpi)
        {
          MPI_Recv(
            i ? &rgt_x0 : &lft_x1,
            1,                              // no of values
            detail::get_mpi_type<real_t>(), // type
            i ? rgt_rank : lft_rank,                       // src comm
            i,                              // message tag
            detail::MPI_COMM_LIBCLOUD,                  // communicator
            MPI_STATUS_IGNORE
          );
        }
      }
#endif
    }
  };
};
