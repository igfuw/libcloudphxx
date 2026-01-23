// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/remove.h>
#include <set>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // remove SDs with n=0
    // loop over all variables, because temporary array of size n is needed,
    // so when done in a single call to remove_if with a large tuple argument,
    // it causes huge memory spikes
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_remove_n0()
    {
      namespace arg = thrust::placeholders;

      // remove chemical stuff
      if(opts_init.chem_switch)
      {
        for (int i = chem_all-1; i >= 0; --i)
        {
          typename thrust_device::vector<real_t>::iterator new_last = thrust::remove_if(
            chem_bgn[i],
            chem_end[i],
            n.begin(),
            arg::_1 == 0
          );
          
          thrust_device::vector<real_t> &vec(
            i < chem_rhs_beg 
              ? chem_ante_rhs
              : i < chem_rhs_fin
                ? chem_rhs
                : chem_post_rhs
          );
 
          vec.erase(new_last, chem_end[i]);
        }
      }

      // remove from real_t vectors
      for(auto vec: distmem_real_vctrs)
      { 
        thrust::remove_if(
          vec.first->begin(),
          vec.first->begin() + n_part,
          n.begin(),
          arg::_1 == 0
        );
      }

      // remove from n and set new n_part
      auto new_end = thrust::remove_if(
        n.begin(),
        n.begin() + n_part,
        arg::_1 == 0
      );
      n_part = new_end - n.begin();

      // resize vectors
      hskpng_resize_npart();

      // resize chem vectors and update chem iterators
      if(opts_init.chem_switch)
        init_chem();

      // NOTE: hskpng_ijk needs to be called aftewards, as i,j,k,ijk were not removed, but were resized - they are undefined!
      // TODO: call it here? bu remove_n0 is called multiple times in rcyc(), so more optimal to do it manually afterwards...
    }
  };  
};
