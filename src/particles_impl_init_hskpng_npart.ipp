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
      if (opts_init.nx != 0) i.resize(opts_init.n_sd_max); //
      if (opts_init.ny != 0) j.resize(opts_init.n_sd_max); //  > TODO: are they needed at all?
      if (opts_init.nz != 0) k.resize(opts_init.n_sd_max); //
      ijk.resize(opts_init.n_sd_max);
      if (n_dims == 0) thrust::fill(ijk.begin(), ijk.end(), 0);

      if (opts_init.nx != 0) x.resize(opts_init.n_sd_max); 
      if (opts_init.ny != 0) y.resize(opts_init.n_sd_max); 
      if (opts_init.nz != 0) z.resize(opts_init.n_sd_max); 

      vt.resize(opts_init.n_sd_max);
      thrust::fill(vt.begin(), vt.end(), 0); // so that it may be safely used in condensation before first update

      sorted_id.resize(opts_init.n_sd_max);
      sorted_ijk.resize(opts_init.n_sd_max);
      
      tmp_device_real_part.resize(opts_init.n_sd_max);
      tmp_device_n_part.resize(opts_init.n_sd_max);

      rd3.resize(opts_init.n_sd_max);
      rw2.resize(opts_init.n_sd_max);
      n.resize(opts_init.n_sd_max);
      kpa.resize(opts_init.n_sd_max);

      sd_stat.resize(opts_init.n_sd_max);
      thrust::fill(sd_stat.begin(), sd_stat.end(), detail::inactive);
    }
  };
};
