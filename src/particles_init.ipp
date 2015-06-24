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
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::init(
      const arrinfo_t<real_t> th,
      const arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> courant_x, // might be NULL
      const arrinfo_t<real_t> courant_y, // might be NULL
      const arrinfo_t<real_t> courant_z  // might be NULL
    )
    {
      // sanity checks
      assert(!th.is_null());
      assert(!rv.is_null());
      assert(!rhod.is_null());

      // initialising Eulerian-Lagrandian coupling
      pimpl->init_sync();
      pimpl->init_e2l(th,   &pimpl->th);
      pimpl->init_e2l(rv,   &pimpl->rv);
      pimpl->init_e2l(rhod, &pimpl->rhod);

      if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0);
      if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0);
      if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1);

      // feeding in Eulerian fields
      pimpl->sync(th,   pimpl->th);
      pimpl->sync(rv,   pimpl->rv);
      pimpl->sync(rhod, pimpl->rhod);
      if (!courant_x.is_null()) pimpl->sync(courant_x, pimpl->courant_x);
      if (!courant_y.is_null()) pimpl->sync(courant_y, pimpl->courant_y);
      if (!courant_z.is_null()) pimpl->sync(courant_z, pimpl->courant_z);

      // initialising housekeeping data of the size ncell
      pimpl->init_hskpng_ncell(); 

      // initialising helper data for advection (Arakawa-C grid neighbours' indices)
      // done before init_xyz, cause it uses dv initialized here
      pimpl->init_grid();

      // initialising particle positions
      pimpl->init_xyz();

      // initialising more housekeeping data (incl. ijk)
      pimpl->hskpng_Tpr(); 
      pimpl->hskpng_ijk(); 

      // initialising dry radii (needs positions, ijk and rhod)
      assert(pimpl->opts_init.dry_distros.size() == 1); // TODO: handle multiple spectra/kappas
      pimpl->init_dry(
        pimpl->opts_init.dry_distros.begin()->first,
        pimpl->opts_init.dry_distros.begin()->second 
      ); // TODO: document that n_of_lnrd_stp is expected!

      // initialising wet radii
      pimpl->init_wet();
      debug::print(pimpl->n);


      // initialising chem stuff
      if(pimpl->opts_init.chem_switch) pimpl->init_chem();

      pimpl->init_sstp();

      //initialising collision kernel
      pimpl->init_kernel();
    }
  };
};
