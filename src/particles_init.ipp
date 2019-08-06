// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */
#include<thrust/extrema.h>

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
      const arrinfo_t<real_t> p,         // pressure profile [in Pascals], needed if pressure perturbations are neglected in condensation (e.g. anelastic model)
                                         // defaults to NULL-NULL pair (variable pressure)
      const arrinfo_t<real_t> courant_x, // might be NULL
      const arrinfo_t<real_t> courant_y, // might be NULL
      const arrinfo_t<real_t> courant_z, // might be NULL
      const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem
    )
    {

      pimpl->init_sanity_check(th, rv, rhod, p, courant_x, courant_y, courant_z, ambient_chem);

      // is a constant pressure profile used?
      pimpl->const_p = !p.is_null();
      // if pressure comes from a profile, sstp_tmp_p also needs to be copied between distributed memories
      if(pimpl->const_p && pimpl->opts_init.sstp_cond > 1 && pimpl->opts_init.exact_sstp_cond)
        pimpl->distmem_real_vctrs.insert(&pimpl->sstp_tmp_p);

      // initialising Eulerian-Lagrangian coupling
      pimpl->init_sync();  // also, init of ambient_chem vectors
      pimpl->init_e2l(th,   &pimpl->th);
      pimpl->init_e2l(rv,   &pimpl->rv);
      pimpl->init_e2l(rhod, &pimpl->rhod);
      if(pimpl->const_p)
        pimpl->init_e2l(p, &pimpl->p);

#if !defined(__NVCC__)
      using std::max;
#endif
      if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0, - pimpl->halo_x);
      if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0, pimpl->n_x_bfr * pimpl->opts_init.nz - pimpl->halo_y);
      if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1, pimpl->n_x_bfr * max(1, pimpl->opts_init.ny) - pimpl->halo_z);
      /*
      if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0);
      if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0, pimpl->n_x_bfr * pimpl->opts_init.nz );
      if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1, pimpl->n_x_bfr * max(1, pimpl->opts_init.ny) );
      */

      if (pimpl->opts_init.chem_switch)
	for (int i = 0; i < chem_gas_n; ++i)
	  pimpl->init_e2l(ambient_chem.at((chem_species_t)i), &pimpl->ambient_chem[(chem_species_t)i]);

      // feeding in Eulerian fields
      pimpl->sync(th,   pimpl->th);
      pimpl->sync(rv,   pimpl->rv);
      pimpl->sync(rhod, pimpl->rhod);
      pimpl->sync(p,   pimpl->p);

      if (!courant_x.is_null()) pimpl->sync(courant_x, pimpl->courant_x, pimpl->halo_x, pimpl->halo_x);
      if (!courant_y.is_null()) pimpl->sync(courant_y, pimpl->courant_y, pimpl->halo_y, pimpl->halo_y);
      if (!courant_z.is_null()) pimpl->sync(courant_z, pimpl->courant_z, pimpl->halo_z, pimpl->halo_z);

      if (pimpl->opts_init.chem_switch)
	for (int i = 0; i < chem_gas_n; ++i)
	  pimpl->sync(
            ambient_chem.at((chem_species_t)i), 
            pimpl->ambient_chem[(chem_species_t)i]
          );

      // initialising housekeeping data of the size ncell
      pimpl->init_hskpng_ncell(); 

      // initialising helper data for advection (Arakawa-C grid neighbours' indices)
      // and cell volumes
      pimpl->init_grid();

      // initialising Tpr
      pimpl->hskpng_Tpr(); 

      // --------  init super-droplets --------
      // reserve memory for data of the size of the max number of SDs
      pimpl->init_hskpng_npart(); 

      // initial parameters (from dry distribution or dry radius-concentration pairs)
      if(pimpl->opts_init.dry_distros.size() > 0)
        pimpl->init_SD_with_distros();
      if(pimpl->opts_init.dry_sizes.size() > 0)
        pimpl->init_SD_with_sizes();

      // --------  other inits  --------
      //initialising collision kernel
      if(pimpl->opts_init.coal_switch) pimpl->init_kernel();

      //initialising vterm
      if(pimpl->opts_init.coal_switch || pimpl->opts_init.sedi_switch) pimpl->init_vterm();

      // initialising neighbouring node's domain sizes
      //if(!pimpl->opts_init.dev_count)
      pimpl->xchng_domains();

      // init count_num and count_ijk
      pimpl->hskpng_count();
    }
  };
};
