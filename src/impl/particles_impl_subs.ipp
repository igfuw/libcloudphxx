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
    void particles_t<real_t, device>::impl::subs()
    {   
      namespace arg = thrust::placeholders;
 
      // settling due to sedimentation + large-scale subsidence
      thrust::transform(
        z.begin(), z.end(),                    // position
        thrust::make_permutation_iterator(w_LS.begin(), k.begin()),     // large-scale subsidence velocity
        z.begin(),                         // output
        arg::_1 - opts_init.dt * arg::_2   // Euler scheme (assuming w_LS positive for downoward motion!)
      );
    }
  };  
};
