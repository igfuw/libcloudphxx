// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    // init
    template <typename real_t, int device>
    void particles<real_t, device>::init()
    {
std::cerr << "\n\n INIT \n\n";

      assert(pimpl->opts.dry_distros.size() == 1); // TODO
      pimpl->init_dry(pimpl->opts.dry_distros.begin()->second); // TODO: document that n_of_lnrd is expected!

      pimpl->init_xyz();
      pimpl->init_Tpr(); // only alloc here?
      pimpl->hskpng(); // could be part of sync_in?
      //pimpl->sync_in();
      pimpl->init_wet();
    }
  };
};
