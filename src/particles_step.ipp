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

      if (pimpl->opts.adve) 
      {
        // advection
        ;
        // periodic boundaries
        ;
      }

      if (pimpl->opts.sedi)
      {
        // sedimentation
        ;
        // recycling
        if (pimpl->opts.rcyc) ;
      }

      // updating particle->cell look-up table
      if (pimpl->opts.adve || pimpl->opts.sedi) 
        pimpl->hskpng_ijk();

      // syncing in Eulerian fields
      pimpl->sync(rhod_th, pimpl->rhod_th);
      pimpl->sync(rhod_rv, pimpl->rhod_rv);
      pimpl->sync(rhod,    pimpl->rhod);

      // updating Tpr look-up table
      pimpl->hskpng_Tpr();

      // condensation/evaporation
      if (pimpl->opts.cond) ;

      // chemistry
      if (pimpl->opts.chem) ;

      // coalescence
      if (pimpl->opts.coal) ;

      // syncing out
      pimpl->sync(pimpl->rhod_th, rhod_th);
      pimpl->sync(pimpl->rhod_rv, rhod_rv);
      pimpl->sync(pimpl->rhod,    rhod);
    }
  };
};
