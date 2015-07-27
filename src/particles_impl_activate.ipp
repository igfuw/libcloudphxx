namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::activate_SDs()
    {
      namespace arg = thrust::placeholders;
      thrust::transform_if( 
        sd_stat.begin(), sd_stat.end(), // input 
        sd_stat.begin(),                // output
        detail::activate<real_t>(),     // operation (make it active)
        detail::to_init()               // condition
      );  
    }  
  };
};
