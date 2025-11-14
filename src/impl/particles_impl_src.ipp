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
    void particles_t<real_t, device>::impl::src(const src_dry_distros_t<real_t> &sdd, const src_dry_sizes_t<real_t> &sds)
    {   
  //    ante_adding_SD();

      if(!sdd.empty())
        src_dry_distros(sdd);

      if(!sds.empty())
        src_dry_sizes(sds);

//      post_adding_SD();
    }
  };  
};
