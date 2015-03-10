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
          k_golovin.resize(1);
          p_kernel =  (&(k_golovin[0])).get();
          break;

        case(geometric):
          k_geometric.resize(1);
          p_kernel =  (&(k_geometric[0])).get();
          break;
          
        default:
          throw std::runtime_error("please supply a type of collision kernel to use");
      }
    }
  }
}
