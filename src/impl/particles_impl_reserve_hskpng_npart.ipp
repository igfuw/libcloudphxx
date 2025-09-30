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
    void particles_t<real_t, device>::impl::reserve_hskpng_npart()
    {
      // memory allocation
      if (opts_init.nx != 0) i.reserve(opts_init.n_sd_max); //
      if (opts_init.ny != 0) j.reserve(opts_init.n_sd_max); //  > TODO: are they needed at all?
      if (opts_init.nz != 0) k.reserve(opts_init.n_sd_max); //
      ijk.reserve(opts_init.n_sd_max);

      if (opts_init.nx != 0) x.reserve(opts_init.n_sd_max); 
      if (opts_init.ny != 0) y.reserve(opts_init.n_sd_max); 
      if (opts_init.nz != 0) z.reserve(opts_init.n_sd_max); 

      if(opts_init.turb_adve_switch)
      {
        if (opts_init.nx != 0) up.reserve(opts_init.n_sd_max);
        if (opts_init.ny != 0) vp.reserve(opts_init.n_sd_max); 
        if (opts_init.nz != 0) wp.reserve(opts_init.n_sd_max); 
      }

      if(opts_init.turb_cond_switch)
      {
        wp.reserve(opts_init.n_sd_max);
        ssp.reserve(opts_init.n_sd_max);
        dot_ssp.reserve(opts_init.n_sd_max);
      }

      if(opts_init.ice_switch)
      {
        ice.reserve(opts_init.n_sd_max);
        rd2_insol.reserve(opts_init.n_sd_max);
        T_freeze.reserve(opts_init.n_sd_max);
        ice_a.reserve(opts_init.n_sd_max);
        ice_c.reserve(opts_init.n_sd_max);
        ice_rho.reserve(opts_init.n_sd_max);
      }

      vt.reserve(opts_init.n_sd_max);

      sorted_id.reserve(opts_init.n_sd_max);
      sorted_ijk.reserve(opts_init.n_sd_max);
      
      tmp_device_real_part.reserve(opts_init.n_sd_max);
      tmp_device_n_part.reserve(opts_init.n_sd_max);
      tmp_device_size_part.reserve(opts_init.n_sd_max);

      rd3.reserve(opts_init.n_sd_max);
      rw2.reserve(opts_init.n_sd_max);
      n.reserve(opts_init.n_sd_max);
      kpa.reserve(opts_init.n_sd_max);

      if(opts_init.diag_incloud_time)
        incloud_time.reserve(opts_init.n_sd_max);

      {
        tmp_device_real_part1.reserve(opts_init.n_sd_max); 
      }
      if((allow_sstp_cond && opts_init.exact_sstp_cond) || n_dims==3 || opts_init.turb_cond_switch)
      {
        tmp_device_real_part2.reserve(opts_init.n_sd_max); 
      }
      if(allow_sstp_cond && opts_init.exact_sstp_cond)
      {
        tmp_device_real_part3.reserve(opts_init.n_sd_max); 
        tmp_device_real_part4.reserve(opts_init.n_sd_max);  
        sstp_tmp_rv.reserve(opts_init.n_sd_max);
        sstp_tmp_th.reserve(opts_init.n_sd_max);
        sstp_tmp_rh.reserve(opts_init.n_sd_max);
        if(const_p) // in const_p pressure is not diagnostic (it's constant) - in per-particle sub-stepping it has to be substepped and we need two vectors to do that
        {
          sstp_tmp_p.reserve(opts_init.n_sd_max);
          tmp_device_real_part5.reserve(opts_init.n_sd_max);  
        }
      }
      // reserve memory for in/out buffers
      // for courant_x = 0.1 and n_sd_max, overkill?
      // done using resize, because _bfr.end() is never used and we want to assert that buffer is large enough using the .size() function
      if(distmem())
      {
        const int no_of_n_vctrs_copied(distmem_n_vctrs.size());
        const int no_of_real_vctrs_copied(distmem_real_vctrs.size());

        in_n_bfr.resize(no_of_n_vctrs_copied * opts_init.n_sd_max / opts_init.nx / config.bfr_fraction);     // for n
        out_n_bfr.resize(no_of_n_vctrs_copied * opts_init.n_sd_max / opts_init.nx / config.bfr_fraction);

        in_real_bfr.resize(no_of_real_vctrs_copied * opts_init.n_sd_max / opts_init.nx / config.bfr_fraction);     // for rd3 rw2 kpa vt x y z  sstp_tmp_th/rv/rh/p, etc.
        out_real_bfr.resize(no_of_real_vctrs_copied * opts_init.n_sd_max / opts_init.nx / config.bfr_fraction);
      }

    // -------- inits done here before resize and reserve were separated. Left for debugging reasons. -----------

//      if (n_dims == 0) thrust::fill(ijk.begin(), ijk.end(), 0);
//
//      // using resize instead of reserve is ok, because hskpng_resize is called right after this init
//      if(opts_init.turb_adve_switch)
//      {
//        if (opts_init.nx != 0) up.resize(opts_init.n_sd_max, 0.); // init with no perturbation 
//        if (opts_init.ny != 0) vp.resize(opts_init.n_sd_max, 0.); 
//        if (opts_init.nz != 0) wp.resize(opts_init.n_sd_max, 0.); 
//      }
//      if(opts_init.turb_cond_switch)
//      {
//        wp.resize(opts_init.n_sd_max, 0.);
//        ssp.resize(opts_init.n_sd_max, 0.);
//        dot_ssp.resize(opts_init.n_sd_max, 0.);
//      }
//      vt.resize(opts_init.n_sd_max, 0.); // so that it may be safely used in condensation before first update
//
//      if(allow_sstp_cond && opts_init.exact_sstp_cond)
//      {
//        sstp_tmp_rv.resize(opts_init.n_sd_max);
//        sstp_tmp_th.resize(opts_init.n_sd_max);
//        sstp_tmp_rh.resize(opts_init.n_sd_max);
//        if(const_p) // in const_p pressure is not diagnostic (it's constant) - in per-particle sub-stepping it has to be substepped and we need two vectors to do that
//        {
//          sstp_tmp_p.resize(opts_init.n_sd_max);
//          tmp_device_real_part5.reserve(opts_init.n_sd_max);  
//        }
//      }
    }
  };
};
