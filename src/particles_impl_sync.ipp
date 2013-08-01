// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include "detail/functors_host.hpp"
#include <thrust/iterator/discard_iterator.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, int device>
    void particles<real_t, device>::impl::sync(
      const arrinfo_t<real_t> &from,
      thrust_device::vector<real_t> &to
    )
    {   
      if (from.is_null()) return;

      thrust::transform(
        l2e[&to].begin(), l2e[&to].end(),
#if defined(__NVCC__) // TODO: better condition (same addressing space)
        tmp_host_real_cell.begin(), 
#else
        to.begin(),
#endif
        detail::c_arr_get<real_t>(from.dataZero)
      );

#if defined(__NVCC__)
      thrust::copy(tmp_host_real_cell.begin(), tmp_host_real_cell.end(), to.begin());
#endif
    }   

    template <typename real_t, int device>
    void particles<real_t, device>::impl::sync(
      const thrust_device::vector<real_t> &from,
      arrinfo_t<real_t> &to
    )
    {   
      if (to.is_null()) return;

#if defined(__NVCC__)
      thrust::copy(from.begin(), from.end(), tmp_host_real_cell.begin());
#endif

      thrust::transform(
        l2e[&from].begin(), l2e[&from].end(), 
#if defined(__NVCC__) // TODO: better condition (same addressing space)
        tmp_host_real_cell.begin(), 
#else
        from.begin(),
#endif
        thrust::make_discard_iterator(),
        detail::c_arr_set<real_t>(to.dataZero)
      );
    }   
  };  
};
