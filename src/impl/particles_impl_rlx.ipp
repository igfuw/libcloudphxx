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
    void particles_t<real_t, device>::impl::rlx()
    {   
      ante_adding_SD();

      if(!opts_init.rlx_dry_distros.empty())
        rlx_dry_distros();
 
      post_adding_SD();
    }
  };  
};
