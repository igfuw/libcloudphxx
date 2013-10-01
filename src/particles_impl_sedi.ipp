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
    void particles_t<real_t, device>::impl::sedi()
    {   
      using namespace thrust::placeholders;
 
      thrust::transform(
        z.begin(), z.end(), // input - 1st arg
        vt.begin(),         // input - 2nd arg
        z.begin(),          // output
        _1 - opts_init.dt * _2   // Euler scheme (assuming vt positive!)
      );
    }
  };  
};
