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
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_sync()
    {
      // memory allocation for scalar fields
//      rhod.resize(n_cell.get());
//      p.resize(n_cell);
//      th.resize(n_cell);
//      rv.resize(n_cell);
      if(opts_init.chem_switch)
        for (int i = 0; i < chem_gas_n; ++i)
          ambient_chem[(chem_species_t)i].resize(n_cell.get());
      if(opts_init.turb_cond_switch || opts_init.turb_adve_switch || opts_init.turb_coal_switch)
        diss_rate.resize(n_cell.get());

      // memory allocation for vector fields (Arakawa-C grid)
      // halo in the x dimensions (which could be distmem boundary)
      switch (n_dims)
      {
        case 3: 
          courant_x.resize((opts_init.nx + 2 * halo_size + 1) * opts_init.ny * opts_init.nz);
          courant_y.resize((opts_init.nx + 2 * halo_size) * (opts_init.ny + 1) * opts_init.nz);
          courant_z.resize((opts_init.nx + 2 * halo_size) * opts_init.ny * (opts_init.nz + 1));
          break;
        case 2: 
          courant_x.resize((opts_init.nx + 2 * halo_size + 1) * opts_init.nz);
          courant_z.resize((opts_init.nx + 2 * halo_size) * (opts_init.nz + 1));
          break;
        case 1:
          courant_x.resize(opts_init.nx + 2 * halo_size + 1);
          break;
        case 0: break;
        default: assert(false); 
      }
    }
  };
};
