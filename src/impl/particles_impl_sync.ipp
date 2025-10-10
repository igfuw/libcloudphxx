// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/iterator/discard_iterator.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sync(
      const arrinfo_t<real_t> &from,
      thrust_device::vector<real_t> &to
    )
    {   
      if (from.is_null()) return;

      assert(to.size() >= l2e[&to].size());

      auto host_consec_g = tmp_host_real_grid.get_guard();
      thrust::host_vector<real_t> &host_consec = host_consec_g.get();

      thrust::transform(
        l2e[&to].begin(), l2e[&to].end(),
#if defined(__NVCC__) // TODO: better condition (same addressing space)
        host_consec.begin(), 
#else
        to.begin(),
#endif
        detail::c_arr_get<real_t>(from.data)
      );

#if defined(__NVCC__)
      thrust::copy(host_consec.begin(), host_consec.begin() + l2e[&to].size(), to.begin());
#endif
    }   

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sync(
      const thrust_device::vector<real_t> &from,
      arrinfo_t<real_t> &to
    )
    {   
      if (to.is_null()) return;

      auto host_consec_g = tmp_host_real_grid.get_guard();
      thrust::host_vector<real_t> &host_consec = host_consec_g.get();

#if defined(__NVCC__)
      assert(from.size() <= host_consec.size());
      thrust::copy(from.begin(), from.end(), host_consec.begin());
#endif

      thrust::transform(
        l2e[&from].begin(), l2e[&from].end(), 
#if defined(__NVCC__) // TODO: better condition (same addressing space)
        host_consec.begin(), 
#else
        from.begin(),
#endif
        thrust::make_discard_iterator(),
        detail::c_arr_set<real_t>(to.data)
      );
    }   
  };  
};
