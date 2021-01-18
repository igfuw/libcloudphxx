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

      thrust::transform(
        l2e[&to].begin(), l2e[&to].end(),
#if defined(__NVCC__) // TODO: better condition (same addressing space)
        tmp_host_real_grid.begin(), 
#else
        to.begin(),
#endif
        detail::c_arr_get<real_t>(from.data)
      );

#if defined(__NVCC__)
      thrust::copy(tmp_host_real_grid.begin(), tmp_host_real_grid.begin() + l2e[&to].size(), to.begin());
#endif
    }   

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sync(
      const thrust_device::vector<real_t> &from,
      arrinfo_t<real_t> &to
    )
    {   
      if (to.is_null()) return;

#if defined(__NVCC__)
      assert(from.size() <= tmp_host_real_grid.size());
      thrust::copy(from.begin(), from.end(), tmp_host_real_grid.begin());
#endif

      thrust::transform(
        l2e[&from].begin(), l2e[&from].end(), 
#if defined(__NVCC__) // TODO: better condition (same addressing space)
        tmp_host_real_grid.begin(), 
#else
        from.begin(),
#endif
        thrust::make_discard_iterator(),
        detail::c_arr_set<real_t>(to.data)
      );
    }   
  };  
};
