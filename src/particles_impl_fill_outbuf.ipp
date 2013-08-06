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
    template <typename real_t, int device>
    void particles<real_t, device>::impl::fill_outbuf(int n)
    {
      thrust::fill(tmp_host_real_cell.begin(), tmp_host_real_cell.end(), 0);

#if defined(__NVCC__)
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

std::cerr << n << " vs. " << outbuf_default_arg << std::endl;

      if (n != outbuf_default_arg)
      { // copying mom (vector of real_t)
	thrust::copy(
	  count_mom[n].begin(),               // input - begin
	  count_mom[n].begin() + count_n,     // input - end
	  thrust::make_permutation_iterator(  // output
	    tmp_host_real_cell.begin(),         // data
	    pi.begin()                          // permutation
	  )
	);
      } 
      else 
      { // copying num (vector of n_t)
	thrust::copy(
	  count_num.begin(),               // input - begin
	  count_num.begin() + count_n,     // input - end
	  thrust::make_permutation_iterator(  // output
	    tmp_host_real_cell.begin(),         // data
	    pi.begin()                          // permutation
	  )
	);
      }
    }
  };
};
