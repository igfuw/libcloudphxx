// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief timestepping routine for super droplets
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    // init
    template <typename real_t, int device>
    void particles<real_t, device>::step(
      real_t *rhod_th,
      real_t *rhod_rv,
      real_t *rhod // defaults to NULL (e.g. kinematic model)
    )
    {
std::cerr << "\n\n STEP \n\n";
      assert(rhod_th != NULL);
      assert(rhod_rv != NULL);

      // syncing in
      pimpl->sync(rhod_th, pimpl->rhod_th);
      pimpl->sync(rhod_rv, pimpl->rhod_rv);
      pimpl->sync(rhod,    pimpl->rhod);

      //
      pimpl->hskpng_ijk();

      // syncing out
      pimpl->sync(pimpl->rhod_th, rhod_th);
      pimpl->sync(pimpl->rhod_rv, rhod_rv);
      pimpl->sync(pimpl->rhod,    rhod);
    }
  };
};
