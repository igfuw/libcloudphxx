namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_kernel()
    {
      switch(opts_init.kernel)
      {
        case(golovin):
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

        case(geometric):
          // init kernel parameters vector
          if(n_kernel_params > 1)
          {
            throw std::runtime_error("Geometric kernel accepts up to one parameter.");
          }
          kernel_parameters.resize(1);
          if(n_kernel_params == 1)
            thrust::copy(opts_init.kernel_parameters.begin(), opts_init.kernel_parameters.end(), kernel_parameters.begin());
          else kernel_parameters[0] = 1.; //default multiplier = 1

          // init kernel
          k_geometric.resize(1, kernel_geometric<real_t, n_t> (kernel_parameters.data()));
          p_kernel = (&(k_geometric[0])).get();
          break;
          
        default:
          throw std::runtime_error("please supply a type of collision kernel to use");
      }
    }
  }
}
