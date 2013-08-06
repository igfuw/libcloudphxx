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
    // records super-droplet concentration per grid cell
    template <typename real_t, int device>
    void particles<real_t, device>::diag_sd_conc()
    {
      pimpl->hskpng_count(); // common code with coalescence, hence separated into a method
    }

    // 
    template <typename real_t, int device>
    void particles<real_t, device>::diag_rw_moms(int)
    {
    }

    //
    template <typename real_t, int device>
    void particles<real_t, device>::diag_rd_moms(int n)
    {
      // locating the requested range
      typename opts_t<real_t>::outmom_t::const_iterator it = pimpl->opts.out_dry.begin();
      while (n-- > 0) it++; 

      pimpl->hskpng_rd_moms(); 
    }
  };
};
