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
      thrust::copy(
        count_ijk.begin(), count_ijk.end(), // from
        tmp_host_size_cell.begin()
      );
#endif

#if !defined(__NVCC__)
      thrust::device_vector<thrust_size_t> &pi(count_ijk);
#else
      thrust::host_vector<thrust_size_t> &pi(tmp_host_size_cell);
#endif

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
