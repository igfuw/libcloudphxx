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
//          k_geometric = thrust_device::vector<kernel_geometric<real_t, n_t> >(1);
          k_geometric.resize(1);
          p_kernel = thrust::constant_iterator<kernel_base<real_t, n_t> *> ( (&(k_geometric[0])).get() );
//          k_golovin.clear();
//          k_golovin.shrink_to_fit();
          break;

        case(golovin):
//          k_golovin = thrust_device::vector<kernel_golovin<real_t, n_t> >(1);
          k_golovin.resize(1);
          p_kernel = thrust::constant_iterator<kernel_base<real_t, n_t> *> ( (&(k_golovin[0])).get() );
//          k_geometric.clear();
//          k_geometric.shrink_to_fit();
          break;

        default:
//          k_geometric = thrust_device::vector<kernel_geometric<real_t, n_t> >(1);
          k_geometric.resize(1);
          p_kernel = thrust::constant_iterator<kernel_base<real_t, n_t> *> ( (&(k_geometric[0])).get() );
//          k_golovin.clear();
//          k_golovin.shrink_to_fit();
          break;
      }
    }
  }
}
