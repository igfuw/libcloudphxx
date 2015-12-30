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
    void particles_t<real_t, device>::impl::init_hskpng_npart()
    {
      // memory allocation
      if (opts_init.nx != 0) i.reserve(opts_init.n_sd_max); //
      if (opts_init.ny != 0) j.reserve(opts_init.n_sd_max); //  > TODO: are they needed at all?
      if (opts_init.nz != 0) k.reserve(opts_init.n_sd_max); //
      ijk.reserve(opts_init.n_sd_max);
      if (n_dims == 0) thrust::fill(ijk.begin(), ijk.end(), 0);

      if (opts_init.nx != 0) x.reserve(opts_init.n_sd_max); 
      if (opts_init.ny != 0) y.reserve(opts_init.n_sd_max); 
      if (opts_init.nz != 0) z.reserve(opts_init.n_sd_max); 

      vt.reserve(opts_init.n_sd_max);
      thrust::fill(vt.begin(), vt.end(), 0); // so that it may be safely used in condensation before first update

      sorted_id.reserve(opts_init.n_sd_max);
      sorted_ijk.reserve(opts_init.n_sd_max);
      
      tmp_device_real_part.reserve(opts_init.n_sd_max);
      tmp_device_n_part.reserve(opts_init.n_sd_max);
      tmp_device_size_part.reserve(opts_init.n_sd_max);

      rd3.reserve(opts_init.n_sd_max);
      rw2.reserve(opts_init.n_sd_max);
      n.reserve(opts_init.n_sd_max);
      kpa.reserve(opts_init.n_sd_max);

      // reserve memory for in/out buffers
      // for courant_x = 0.1 and n_sd_max
      // overkill?
      if(dist_mem)
      {
        in_n_bfr.resize(opts_init.n_sd_max / opts_init.nx / 10);     // for n
        out_n_bfr.resize(opts_init.n_sd_max / opts_init.nx / 10);
        in_real_bfr.resize(7 * opts_init.n_sd_max / opts_init.nx / 10);     // for rd3 rw2 kpa vt x y z
        out_real_bfr.resize(7 * opts_init.n_sd_max / opts_init.nx / 10);
      }
      if(opts_init.chem_switch)
      {
        tmp_device_real_part_chem.reserve(opts_init.n_sd_max); 
        V_old.reserve(opts_init.n_sd_max);
        tmp_device_real_part_HNO3.reserve(opts_init.n_sd_max); // TODO: only in chemistry, but can we do it without?
        tmp_device_real_part_NH3.reserve(opts_init.n_sd_max);  // TODO: only in chemistry, but can we do it without?
        tmp_device_real_part_CO2.reserve(opts_init.n_sd_max);  // TODO: only in chemistry, but can we do it without?
        tmp_device_real_part_SO2.reserve(opts_init.n_sd_max);  // TODO: only in chemistry, but can we do it without?
      }
    }
  };
};
