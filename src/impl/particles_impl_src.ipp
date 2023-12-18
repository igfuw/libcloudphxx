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
    void particles_t<real_t, device>::impl::src(const real_t &dt)
    {   
  //    ante_adding_SD();

      if(!opts.src_dry_distros.empty())
        src_dry_distros(dt);

      if(!opts.src_dry_sizes.empty())
        src_dry_sizes(dt);

//      post_adding_SD();
    }
  };  
};
