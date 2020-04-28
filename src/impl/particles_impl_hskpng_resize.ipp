namespace libcloudphxx
{
  namespace lgrngn
  {
    // resize vectors to n_part
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_resize_npart()
    {
      if(n_part > opts_init.n_sd_max) throw std::runtime_error(detail::formatter() << "n_sd_max (" << opts_init.n_sd_max << ") < n_part (" << n_part << ")");
      {
        thrust_device::vector<real_t> *vec_arr[] = {&rw2, &rd3, &kpa, &tmp_device_real_part, &delta_revp20, &delta_revp25, &delta_revp32, &delta_accr20, &delta_accr25, &delta_acnv20, &delta_acnv25, &delta_accr32, &delta_acnv32};
        for(auto vec : vec_arr)
        {
          vec->resize(n_part);
        }
      }
      {
        thrust_device::vector<thrust_size_t> *vec_arr[] = {&ijk, &sorted_id, &sorted_ijk};
        for(auto vec : vec_arr)
        {
          vec->resize(n_part);
        }
      }
      n.resize(n_part);
      tmp_device_n_part.resize(n_part);
      tmp_device_size_part.resize(n_part);

      vt.resize(n_part, detail::invalid);

      if (opts_init.nx != 0) i.resize(n_part); 
      if (opts_init.ny != 0) j.resize(n_part); 
      if (opts_init.nz != 0) k.resize(n_part); 

      if (opts_init.nx != 0) x.resize(n_part); 
      if (opts_init.ny != 0) y.resize(n_part); 
      if (opts_init.nz != 0) z.resize(n_part); 

      if(opts_init.turb_adve_switch) 
      {
        if (opts_init.nx != 0) up.resize(n_part, 0); 
        if (opts_init.ny != 0) vp.resize(n_part, 0); 
        if (opts_init.nz != 0) wp.resize(n_part, 0); 
      }

      if(opts_init.turb_cond_switch)
      {
        wp.resize(n_part, 0);
        ssp.resize(n_part, 0);
        dot_ssp.resize(n_part, 0);
      }

      if(opts_init.chem_switch || opts_init.sstp_cond > 1 || n_dims >= 2)
      {
        tmp_device_real_part1.resize(n_part);
      }
      if((opts_init.sstp_cond>1 && opts_init.exact_sstp_cond) || n_dims==3 || opts_init.turb_cond_switch)
      {
        tmp_device_real_part2.resize(n_part);
      }

      if(opts_init.sstp_cond>1 && opts_init.exact_sstp_cond)
      {
        tmp_device_real_part3.resize(n_part);
        tmp_device_real_part4.resize(n_part);  
        sstp_tmp_rv.resize(n_part);
        sstp_tmp_th.resize(n_part);
        sstp_tmp_rh.resize(n_part);
        if(const_p)
        {
          tmp_device_real_part5.resize(n_part);  
          sstp_tmp_p.resize(n_part);
        }
      }

      if(opts_init.diag_incloud_time)
        incloud_time.resize(n_part);
    }
  };
};

