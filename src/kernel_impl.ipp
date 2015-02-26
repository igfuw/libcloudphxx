namespace libcloudphxx
{
  namespace lgrngn
  {
    using detail::tpl_rw_t_wrap;

    template <typename real_t, typename n_t>
    struct kernel_impl_t
    {
//      const std::vector<real_t> params;  kernel parameters from opts_init
      real_t (*pkernel_function)(const tpl_rw_t_wrap<real_t,n_t>& );// TODO: pass parameters

      kernel_impl_t(real_t (*_pkernel_function)(const tpl_rw_t_wrap<real_t,n_t>&) ) : //TODO: pass kernel parameters
//        params(_params),
        pkernel_function(_pkernel_function)
        {}


      BOOST_GPU_ENABLED
      real_t operator() (const tpl_rw_t_wrap<real_t,n_t> &tpl_wrap) const
      {
        return pkernel_function(tpl_wrap); //TODO: pass parameters
      }
    };

    //kernel factory
    template <typename real_t, typename n_t>
    kernel_impl_t<real_t,n_t> kernel_factory(const opts_init_t<real_t> &opts_init)
    {
      switch (opts_init.kernel)
      {
        //TODO: add parameters check
        case(geometric):
          return kernel_impl_t<real_t,n_t>( geometric_kernel<real_t,n_t>);

        case(golovin):
          return kernel_impl_t<real_t,n_t>( golovin_kernel<real_t,n_t>);

        default:
          return kernel_impl_t<real_t,n_t>( geometric_kernel<real_t,n_t>);
      }
    }
  }
}



