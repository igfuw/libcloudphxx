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
    template <typename real_t, int device>
    void particles<real_t, device>::step(
      arrinfo_t<real_t> rhod_th,
      arrinfo_t<real_t> rhod_rv,
      const arrinfo_t<real_t> courant_x, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> courant_y, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> courant_z, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> rhod       // defaults to NULL-NULL pair (e.g. kinematic or boussinesq model)
    )
    {
std::cerr << "\n\n STEP \n\n";

      // syncing in Eulerian fields (if not null)
      pimpl->sync(rhod_th,   pimpl->rhod_th);
      pimpl->sync(rhod_rv,   pimpl->rhod_rv);
      pimpl->sync(courant_x, pimpl->courant_x);
      pimpl->sync(courant_y, pimpl->courant_y);
      pimpl->sync(courant_z, pimpl->courant_z);
      pimpl->sync(rhod,      pimpl->rhod);

      // changing droplet positions
      if (pimpl->opts.adve) 
      {
        // advection
        pimpl->adve();
        // periodic boundaries
        ; // TODO? (if needed)
      }
      if (pimpl->opts.sedi) ; //pimpl->sedi(); // TODO
      if (pimpl->opts.adve || pimpl->opts.sedi) 
      {
        // recycling out-of-domain particles (due to precipitation or advection errors)
        if (pimpl->opts.rcyc) ; // TODO

        // updating particle->cell look-up table
        pimpl->hskpng_ijk();
      }

      // updating Tpr look-up table
      pimpl->hskpng_Tpr();

      // condensation/evaporation
      if (pimpl->opts.cond) ; // TODO

      // chemistry
      if (pimpl->opts.chem) assert(false && "unimplemented");

      // coalescence
      if (pimpl->opts.coal) ; // TODO

      // syncing out
      //pimpl->sync(pimpl->rhod_th, rhod_th);
      //pimpl->sync(pimpl->rhod_rv, rhod_rv);
    }
  };
};
