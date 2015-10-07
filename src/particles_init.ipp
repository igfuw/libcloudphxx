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
      if (pimpl->init_called) 
        throw std::runtime_error("init() may be called just once");
      pimpl->init_called = true;

      // sanity checks
      assert(!th.is_null());
      assert(!rv.is_null());
      assert(!rhod.is_null());

      // --------  init cell characteristics  --------
      // initialising Eulerian-Lagrandian coupling
      pimpl->init_sync();
      pimpl->init_e2l(th,   &pimpl->th);
      pimpl->init_e2l(rv,   &pimpl->rv);
      pimpl->init_e2l(rhod, &pimpl->rhod);

#if !defined(__NVCC__)
      using std::max;
#endif
      if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0);
      if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0, pimpl->n_x_bfr * pimpl->opts_init.nz);
      if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1, pimpl->n_x_bfr * max(1, pimpl->opts_init.ny));

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

      // initialising Tpr
      pimpl->hskpng_Tpr(); 

      pimpl->init_sstp();

      // --------  init super-droplet characteristics  --------
      // reserve memory for data of the size of the max number of SDs
      pimpl->init_hskpng_npart(); 

      // init number of SDs in cells
      pimpl->init_count_num();

      // update no of particles
      // TODO: move to a separate function
      pimpl->n_part_old = 0;
      pimpl->n_part_to_init = thrust::reduce(pimpl->count_num.begin(), pimpl->count_num.end());
      pimpl->n_part += pimpl->n_part_to_init;
      pimpl->hskpng_resize_npart(); 

      // init ijk vector, also n_part and resize n_part vectors
      pimpl->init_ijk();

      // initialising dry radii (needs ijk and rhod)
      assert(pimpl->opts_init.dry_distros.size() == 1); // TODO: handle multiple spectra/kappas
      // analyze the distribution;
      pimpl->dist_analysis(
        pimpl->opts_init.dry_distros.begin()->second,
        pimpl->opts_init.sd_conc
      );
      pimpl->init_dry();

      // init multiplicities
      pimpl->init_n(
        pimpl->opts_init.dry_distros.begin()->first,
        pimpl->opts_init.dry_distros.begin()->second
      ); // TODO: document that n_of_lnrd_stp is expected!

      // initialising wet radii
      pimpl->init_wet();

      // initialising particle positions
      pimpl->init_xyz();

      // initialising chem stuff
      if(pimpl->opts_init.chem_switch) pimpl->init_chem();

      // --------  other inits  --------
      //initialising collision kernel
      if(pimpl->opts_init.coal_switch) pimpl->init_kernel();
    }
  };
};
