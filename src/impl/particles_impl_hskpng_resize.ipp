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
        thrust_device::vector<real_t> *vec[] = {&tmp_device_real_part};
        for(int i=0; i<4; ++i)
        {
          vec[i]->resize(n_part);
        }
      }
      {
        thrust_device::vector<thrust_size_t> *vec[] = {&ijk, &sorted_id, &sorted_ijk};
        for(int i=0; i<3; ++i)
        {
          vec[i]->resize(n_part);
        }
      }
      n.resize(n_part);
      tmp_device_n_part.resize(n_part);
      tmp_device_size_part.resize(n_part);

      if (opts_init.nx != 0) i.resize(n_part); 
      if (opts_init.ny != 0) j.resize(n_part); 
      if (opts_init.nz != 0) k.resize(n_part); 

      if(opts_init.chem_switch || allow_sstp_cond || n_dims >= 2)
      {
        tmp_device_real_part1.resize(n_part);
      }
      if((allow_sstp_cond && opts_init.exact_sstp_cond) || n_dims==3 || opts_init.turb_cond_switch)
      {
        tmp_device_real_part2.resize(n_part);
      }

      if(allow_sstp_cond && opts_init.exact_sstp_cond)
      {
        tmp_device_real_part3.resize(n_part);
        tmp_device_real_part4.resize(n_part);  
        if(const_p)
        {
          tmp_device_real_part5.resize(n_part);  
        }
      }

      for(auto &pair: distmem_real_vctrs)
      {
        if(pair.second == detail::no_initial_value)
          pair.first->resize(n_part);
        else
          pair.first->resize(n_part, pair.second);
      }
    }
  };
};

