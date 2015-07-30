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
        thrust_device::vector<real_t> *vec[] = {&rw2, &rd3, &kpa, &x, &y, &z, &vt, &tmp_device_real_part, &u01};
        for(int i=0; i<9; ++i)
        {
          vec[i]->resize(n_part);
        }
      }
      {
        thrust_device::vector<thrust_size_t> *vec[] = {&i, &j, &k, &ijk, &sorted_id, &sorted_ijk};
        for(int i=0; i<6; ++i)
        {
          vec[i]->resize(n_part);
        }
      }
      n.resize(n_part);
      un.resize(n_part);
    }
  };
};

