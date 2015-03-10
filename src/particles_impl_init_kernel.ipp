namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_kernel()
    {
      switch(opts_init.kernel)
      {
        case(geometric):
          k_geometric.resize(1);
          p_kernel = thrust::constant_iterator<kernel_base<real_t, n_t> *> ( (&(k_geometric[0])).get() );
          break;

        case(golovin):
          k_golovin.resize(1);
          p_kernel = thrust::constant_iterator<kernel_base<real_t, n_t> *> ( (&(k_golovin[0])).get() );
          break;

        default:
          k_geometric.resize(1);
          p_kernel = thrust::constant_iterator<kernel_base<real_t, n_t> *> ( (&(k_geometric[0])).get() );
          break;
      }
    }
  }
}
