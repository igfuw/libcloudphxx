// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include "detail/functors_host.hpp"
#include <thrust/for_each.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, int device>
    void particles<real_t, device>::impl::sync(
      real_t *from,
      thrust_device::vector<real_t> &to
    )
    {   
      if (from == NULL) return;
      struct copy
      {
        real_t *from; // member field
        copy(real_t *from) : from(from) {} // ctor
        real_t operator()(thrust_size_t ix) { return from[ix]; } // op invoked by transform
      };
      thrust::transform(l2e.begin(), l2e.end(), to.begin(), copy(from));
    }   

    template <typename real_t, int device>
    void particles<real_t, device>::impl::sync(
      const thrust_device::vector<real_t> &from,
      real_t *to
    )
    {   
      if (to == NULL) return;
      struct copy
      {
        // member fields
        const thrust_device::vector<real_t> &from;
        real_t *to; 
        const thrust_device::vector<thrust_size_t> &l2e;

        // ctor
        copy(
          const thrust_device::vector<real_t> &from, 
          real_t *to, 
          const thrust_device::vector<thrust_size_t> &l2e
        ) : from(from), to(to), l2e(l2e) {}

        // operator to be invoked by for_each
        void operator()(thrust_size_t ix)
        {
          to[l2e[ix]] = from[ix];
        } 
      };
      thrust::for_each(zero, zero + n_cell, copy(from, to, l2e));
    }   
  };  
};
