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
          // init kernel parameters vector, TODO: parameters sanity check
          kernel_parameters.resize(1, 2000.); //TODO: zamiast 1500 dac 4/3 * b podane w opts_init

          // init kernel
          k_golovin.resize(1, kernel_golovin<real_t, n_t> (kernel_parameters.data()));
          p_kernel = (&(k_golovin[0])).get();
          break;

        case(geometric):
          // init kernel parameters vector, TODO: parameters sanity check

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
