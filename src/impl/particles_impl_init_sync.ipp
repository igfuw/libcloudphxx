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
      rhod.resize(n_cell);
      th.resize(n_cell);
      rv.resize(n_cell);
      RH.resize(n_cell); 
      for (int i = 0; i < chem_gas_n; ++i)
        ambient_chem[(chem_species_t)i].resize(n_cell);

      // memory allocation for vector fields (Arakawa-C grid)
      // 1-cell halo in the x dimensions (which could be distmem boundary)
      switch (n_dims)
      {
        case 3: 
          courant_x.resize((opts_init.nx + 2 + 1) * opts_init.ny * opts_init.nz);
          courant_y.resize((opts_init.nx + 2) * (opts_init.ny + 1) * opts_init.nz);
          courant_z.resize((opts_init.nx + 2) * opts_init.ny * (opts_init.nz + 1));
          break;
        case 2: 
          courant_x.resize((opts_init.nx + 2 + 1) * opts_init.nz);
          courant_z.resize((opts_init.nx + 2) * (opts_init.nz + 1));
          break;
        case 1:
          courant_x.resize(opts_init.nx + 2 + 1);
          break;
        case 0: break;
        default: assert(false); 
      }
    }
  };
};
