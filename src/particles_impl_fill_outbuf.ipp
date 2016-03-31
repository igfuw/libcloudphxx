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
    void particles_t<real_t, device>::impl::fill_outbuf()
    {
      thrust::fill(tmp_host_real_cell.begin(), tmp_host_real_cell.end(), 0);

#if defined(__NVCC__)
printf("rank %d przed copy raz w fill_outbuf\n", mpi_rank);
      thrust::copy(
        count_ijk.begin(), count_ijk.end(), // from
        tmp_host_size_cell.begin()
      );
#endif

#if !defined(__NVCC__)
      thrust_device::vector<thrust_size_t> &pi(count_ijk);
#else
      thrust::host_vector<thrust_size_t> &pi(tmp_host_size_cell);
#endif

printf("rank %d przed copy dwa w fill_outbuf count_n %d\n", mpi_rank, int(count_n));
      thrust::copy(
	count_mom.begin(),               // input - begin
	count_mom.begin() + count_n,     // input - end
	thrust::make_permutation_iterator(  // output
	  tmp_host_real_cell.begin(),         // data
	  pi.begin()                          // permutation
	)
      );
    }
  };
};
