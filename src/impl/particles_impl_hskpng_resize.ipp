namespace libcloudphxx
{
  namespace lgrngn
  {
    // resize vectors to n_part
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_resize_npart()
    {
      if(n_part > opts_init.n_sd_max) throw std::runtime_error(detail::formatter() << "n_sd_max (" << opts_init.n_sd_max << ") < n_part (" << n_part << ")");

      for(auto &vec: distmem_n_vctrs)
        vec->resize(n_part);

      for(auto &pair: distmem_real_vctrs)
      {
        if(pair.second == detail::no_initial_value)
          pair.first->resize(n_part);
        else
          pair.first->resize(n_part, pair.second);
      }

//      for(auto &vec: resize_real_vctrs)
//        vec->resize(n_part);
      tmp_device_real_part.resize(n_part);
      tmp_device_n_part.resize(n_part);
      // tmp_device_size_part.resize(n_part);
      tmp_host_size_part.resize(n_part);
      tmp_host_real_part.resize(n_part);

      for(auto &vec: resize_size_vctrs)
        vec->resize(n_part);

    }
  };
};

