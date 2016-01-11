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

      sstp_tmp_rv.resize(opts_init.n_sd_max);
      sstp_tmp_th.resize(opts_init.n_sd_max);
      sstp_tmp_rh.resize(opts_init.n_sd_max);

      if(opts_init.chem_switch || (opts_init.sstp_cond>1))
      {
        tmp_device_real_part1.reserve(opts_init.n_sd_max); 
        tmp_device_real_part2.reserve(opts_init.n_sd_max); 
        tmp_device_real_part3.reserve(opts_init.n_sd_max); 
        tmp_device_real_part4.reserve(opts_init.n_sd_max);  
        tmp_device_real_part5.reserve(opts_init.n_sd_max);  
        tmp_device_real_part6.reserve(opts_init.n_sd_max);  
      }

      if(opts_init.chem_switch)
      {
      V_old.reserve(opts_init.n_sd_max);
      }
    }
  };
};
