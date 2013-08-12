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

      // updating Tpr look-up table
      pimpl->hskpng_Tpr(); // note: before sedi (needed by v_term)

      // changing droplet positions
      if (pimpl->opts.adve) 
      {
        // advection (TODO temporarily it also includes periodic boundaries)
        pimpl->adve();
        // boundary condition
        ; // TODO
      }
      if (pimpl->opts.sedi) 
      {
        // advection with terminal velocity
        //pimpl->sedi(); // TODO

        // recycling out-of-domain particles (due to precipitation)
        if (pimpl->opts.rcyc) ; // TODO

        // updating particle->cell look-up table
      }
      if (pimpl->opts.adve || pimpl->opts.sedi)
      {
        pimpl->hskpng_ijk();
      }

      // condensation/evaporation
      if (pimpl->opts.cond) for (int step = 0; step < pimpl->opts.sstp_cond; ++step) 
      { 
        pimpl->cond(pimpl->opts.dt / pimpl->opts.sstp_cond, pimpl->opts.RH_max); 
        pimpl->hskpng_Tpr(); // needed even at last iteration (chem & v_term in sedi)
      }

      // chemistry
      if (pimpl->opts.chem) ; // TODO assert(false && "unimplemented");

      // coalescence
      if (pimpl->opts.coal) ; // TODO

      // syncing out (TODO: if cond || ...)
      pimpl->sync(pimpl->rhod_th, rhod_th);
      pimpl->sync(pimpl->rhod_rv, rhod_rv);
    }
  };
};
