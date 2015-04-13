namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_kernel()
    {
      switch(opts_init.kernel)
      {
        case(kernel_t::golovin):
          // init device kernel parameters vector
          if(n_kernel_params != 1)
          {
            throw std::runtime_error("Golovin kernel accepts exactly one parameter.");
          }
          kernel_parameters.resize(n_kernel_params);
          thrust::copy(opts_init.kernel_parameters.begin(), opts_init.kernel_parameters.end(), kernel_parameters.begin());

          // init kernel
          k_golovin.resize(1, kernel_golovin<real_t, n_t> (kernel_parameters.data()));
          p_kernel = (&(k_golovin[0])).get();
          break;

        case(kernel_t::geometric):
          // init kernel parameters vector
          if(n_kernel_params != 0)
          {
            throw std::runtime_error("Geometric kernel doesn't accept parameters.");
          }

          // init kernel
          k_geometric.resize(1, kernel_geometric<real_t, n_t> (kernel_parameters.data()));
          p_kernel = (&(k_geometric[0])).get();
          break;
          
        default:
          ;
      }
    }
  }
}
