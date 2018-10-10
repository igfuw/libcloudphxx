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

      // using resize instead of reserve is ok, because hskpng_resize is called right after this init
      if(opts_init.turb_adve_switch)
      {
        if (opts_init.nx != 0) up.resize(opts_init.n_sd_max, 0.); // init with no perturbation 
        if (opts_init.ny != 0) vp.resize(opts_init.n_sd_max, 0.); 
        if (opts_init.nz != 0) wp.resize(opts_init.n_sd_max, 0.); 
      }

      if(opts_init.turb_cond_switch)
      {
        wp.resize(opts_init.n_sd_max, 0.);
        ssp.resize(opts_init.n_sd_max, 0.);
        dot_ssp.resize(opts_init.n_sd_max, 0.);
      }

      vt.resize(opts_init.n_sd_max, 0.); // so that it may be safely used in condensation before first update

      sorted_id.reserve(opts_init.n_sd_max);
      sorted_ijk.reserve(opts_init.n_sd_max);
      
      tmp_device_real_part.reserve(opts_init.n_sd_max);
      tmp_device_n_part.reserve(opts_init.n_sd_max);
      tmp_device_size_part.reserve(opts_init.n_sd_max);

      rd3.reserve(opts_init.n_sd_max);
      rw2.reserve(opts_init.n_sd_max);
      n.reserve(opts_init.n_sd_max);
      kpa.reserve(opts_init.n_sd_max);

      if(opts_init.chem_switch || opts_init.sstp_cond > 1 || n_dims >= 2)
      {
        tmp_device_real_part1.reserve(opts_init.n_sd_max); 
      }
      if((opts_init.sstp_cond>1 && opts_init.exact_sstp_cond) || n_dims==3 || opts_init.turb_cond_switch)
      {
        tmp_device_real_part2.reserve(opts_init.n_sd_max); 
      }
      if(opts_init.sstp_cond>1 && opts_init.exact_sstp_cond)
      {
        tmp_device_real_part3.reserve(opts_init.n_sd_max); 
        tmp_device_real_part4.reserve(opts_init.n_sd_max);  
        sstp_tmp_rv.resize(opts_init.n_sd_max);
        sstp_tmp_th.resize(opts_init.n_sd_max);
        sstp_tmp_rh.resize(opts_init.n_sd_max);
        if(const_p) // in const_p pressure is not diagnostic (it's constant) - in per-particle sub-stepping it has to be substepped and we need two vectors to do that
        {
          sstp_tmp_p.resize(opts_init.n_sd_max);
          tmp_device_real_part5.reserve(opts_init.n_sd_max);  
        }
      }
      // reserve memory for in/out buffers
      if(opts_init.dev_count > 1)
      {
        in_n_bfr.resize(opts_init.n_sd_max / opts_init.nx / config.bfr_fraction);     // for n
        out_n_bfr.resize(opts_init.n_sd_max / opts_init.nx / config.bfr_fraction);

        in_real_bfr.resize(11 * opts_init.n_sd_max / opts_init.nx / config.bfr_fraction);     // for rd3 rw2 kpa vt x y z  sstp_tmp_th/rv/rh/p
        out_real_bfr.resize(11 * opts_init.n_sd_max / opts_init.nx / config.bfr_fraction);
      }
    }
  };
};
