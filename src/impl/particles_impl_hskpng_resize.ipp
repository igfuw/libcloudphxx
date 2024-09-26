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
      {
//        std::cerr << "resize distmem_n_vctrs vec address: " << vec << std::endl;
//        std::cerr << "vec pre resize: " << std::endl;
//        debug::print(*vec);
        vec->resize(n_part);
//        std::cerr << "vec post resize: " << std::endl;
//        debug::print(*vec);
      }

      for(auto &pair: distmem_real_vctrs)
      {
        if(pair.second == detail::no_initial_value)
          pair.first->resize(n_part);
        else
          pair.first->resize(n_part, pair.second);
      }

      for(auto &vec: resize_real_vctrs)
        vec->resize(n_part);

      for(auto &vec: resize_size_vctrs)
        vec->resize(n_part);

      // its unsigned int vector, probably only one we will use, hence no resize_t_vctrs helper used
      tmp_device_n_part.resize(n_part);

//    not needed, since ijk_history is already in distmem_n_vctrs
//      for(auto &arr : ijk_history)
//        arr.resize(n_part);
    }
  };
};

