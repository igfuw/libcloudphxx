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
      const arrinfo_t<real_t> courant_z, // might be NULL
      const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem
    )
    {
      if (pimpl->init_called) 
        throw std::runtime_error("init() may be called just once");
      pimpl->init_called = true;

      // sanity checks
      if (th.is_null() || rv.is_null() || rhod.is_null())
        throw std::runtime_error("passing th, rv and rhod is mandatory");

      if (!courant_x.is_null() || !courant_y.is_null() || !courant_z.is_null())
      {
	if (pimpl->n_dims == 0)
	  throw std::runtime_error("Courant numbers passed in 0D setup");

	if (pimpl->n_dims == 1 && (courant_x.is_null() || !courant_y.is_null() || !courant_z.is_null()))
	  throw std::runtime_error("Only X Courant number allowed in 1D setup");

	if (pimpl->n_dims == 2 && (courant_x.is_null() || !courant_y.is_null() || courant_z.is_null()))
	  throw std::runtime_error("Only X and Z Courant numbers allowed in 2D setup");

	if (pimpl->n_dims == 3 && (courant_x.is_null() || courant_y.is_null() || courant_z.is_null()))
	  throw std::runtime_error("All XYZ Courant number components required in 3D setup");
      }

      if (pimpl->opts_init.chem_switch && ambient_chem.size() != chem_gas_n) 
        throw std::runtime_error("chemistry was not switched off and ambient_chem is empty");

      if (!pimpl->opts_init.chem_switch && ambient_chem.size() != 0) 
        throw std::runtime_error("chemistry was switched off and ambient_chem is not empty");

      // initialising chem stuff (at the beginning to have ambient_chem vectors allocated)
      if(pimpl->opts_init.chem_switch) pimpl->init_chem();

      // initialising Eulerian-Lagrangian coupling
      pimpl->init_sync();
      pimpl->init_e2l(th,   &pimpl->th);
      pimpl->init_e2l(rv,   &pimpl->rv);
      pimpl->init_e2l(rhod, &pimpl->rhod);

      if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0);
      if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0);
      if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1);

      if (pimpl->opts_init.chem_switch)
	for (int i = 0; i < chem_gas_n; ++i)
	  pimpl->init_e2l(ambient_chem.at((chem_species_t)i), &pimpl->ambient_chem[(chem_species_t)i]);

      // feeding in Eulerian fields
      pimpl->sync(th,   pimpl->th);
      pimpl->sync(rv,   pimpl->rv);
      pimpl->sync(rhod, pimpl->rhod);

      if (!courant_x.is_null()) pimpl->sync(courant_x, pimpl->courant_x);
      if (!courant_y.is_null()) pimpl->sync(courant_y, pimpl->courant_y);
      if (!courant_z.is_null()) pimpl->sync(courant_z, pimpl->courant_z);

      if (pimpl->opts_init.chem_switch)
	for (int i = 0; i < chem_gas_n; ++i)
	  pimpl->sync(
            ambient_chem.at((chem_species_t)i), 
            pimpl->ambient_chem[(chem_species_t)i]
          );

      // initialising particle positions
      pimpl->init_xyz();

      // initialising helper data for advection (Arakawa-C grid neighbours' indices)
      pimpl->init_grid();

      // initialising housekeeping data (incl. ijk)
      pimpl->init_hskpng(); 
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

      // calculate initail volume (helper for Henry in chem)
      if (pimpl->opts_init.chem_switch){
        pimpl->chem_vol_ante();
        pimpl->chem_vol_post();
      }
	
      pimpl->init_sstp();

      //initialising collision kernel
      if(pimpl->opts_init.coal_switch) pimpl->init_kernel();
    }
  };
};
