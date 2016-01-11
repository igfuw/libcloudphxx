namespace libcloudphxx
{
  namespace lgrngn
  {
    // resize vectors to n_part
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_resize_npart()
    {
      if(n_part > opts_init.n_sd_max) throw std::runtime_error("n_sd_max < n_part");
      {
        thrust_device::vector<real_t> *vec[] = {&rw2, &rd3, &kpa, &vt, &tmp_device_real_part};
        for(int i=0; i<5; ++i)
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

      sstp_tmp_rv.resize(n_part);
      sstp_tmp_th.resize(n_part);
      sstp_tmp_rh.resize(n_part);

      if (opts_init.nx != 0) i.resize(n_part); 
      if (opts_init.ny != 0) j.resize(n_part); 
      if (opts_init.nz != 0) k.resize(n_part); 

      if (opts_init.nx != 0) x.resize(n_part); 
      if (opts_init.ny != 0) y.resize(n_part); 
      if (opts_init.nz != 0) z.resize(n_part); 

      if(opts_init.chem_switch || (opts_init.sstp_cond>1))
      {
        tmp_device_real_part1.resize(n_part);
        tmp_device_real_part2.resize(n_part);
        tmp_device_real_part3.resize(n_part);
        tmp_device_real_part4.resize(n_part);  
        tmp_device_real_part5.resize(n_part); 
        tmp_device_real_part6.resize(n_part);
      }
      if(opts_init.chem_switch)
        V_old.resize(n_part);
    }
  };
};

