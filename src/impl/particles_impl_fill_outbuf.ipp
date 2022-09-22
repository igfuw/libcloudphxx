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
      auto &pi(count_ijk);
#else
      auto &pi(tmp_host_size_cell);
#endif

      thrust::copy(
        count_mom.begin(),               // input - begin
        count_mom.begin() + count_n.get(),     // input - end
        thrust::make_permutation_iterator(  // output
          tmp_host_real_cell.begin(),         // data
          pi.begin()                          // permutation
        )
      );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::fill_outbuf_ref()
    {
      thrust::fill(tmp_host_real_cell.begin_ref(), tmp_host_real_cell.end_ref(), 0);

#if defined(__NVCC__)
      thrust::copy(
        count_ijk.begin_ref(), count_ijk.end_ref(), // from
        tmp_host_size_cell.begin_ref()
      );
#endif

#if !defined(__NVCC__)
      auto &pi(count_ijk);
#else
      auto &pi(tmp_host_size_cell);
#endif

      thrust::copy(
        count_mom.begin_ref(),               // input - begin
        count_mom.begin_ref() + count_n.get_ref(),     // input - end
        thrust::make_permutation_iterator(  // output
          tmp_host_real_cell.begin_ref(),         // data
          pi.begin_ref()                          // permutation
        )
      );
    }
  };
};
