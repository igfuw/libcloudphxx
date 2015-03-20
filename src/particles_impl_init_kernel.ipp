namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_kernel()
    {
      //temporary vector to store kernel efficiencies before they are appended to user input parameters
      std::vector<real_t> tmp_kernel_eff;

      switch(opts_init.kernel)
      {
        case(golovin):
          if(n_user_params != 1)
          {
            throw std::runtime_error("Golovin kernel accepts exactly one parameter.");
          }
          // init device kernel parameters vector
          kernel_parameters.resize(n_user_params);
          thrust::copy(opts_init.kernel_parameters.begin(), opts_init.kernel_parameters.end(), kernel_parameters.begin());

          // init kernel
          k_golovin.resize(1, kernel_golovin<real_t, n_t> (kernel_parameters.data()));
          p_kernel = (&(k_golovin[0])).get();
          break;

        case(geometric):
          if(n_user_params != 0)
          {
            throw std::runtime_error("Geometric kernel doesn't accept parameters.");
          }

          // init kernel
          k_geometric.resize(1, kernel_geometric<real_t, n_t> (kernel_parameters.data()));
          p_kernel = (&(k_geometric[0])).get();
          break;
          
        case(hall_davis_no_waals):
          if(n_user_params != 0)
          {
            throw std::runtime_error("Hall + Davis kernel doesn't accept parameters.");
          }
          //read in kernel efficiencies to a temporary container
          detail::hall_davis_no_waals_efficiencies<real_t> (tmp_kernel_eff);
         
          //reserve device memory for kernel parameters vector
          kernel_parameters.resize(opts_init.kernel_parameters.size() + tmp_kernel_eff.size());

          //copy user-defined parameters to device memory
          thrust::copy(opts_init.kernel_parameters.begin(), opts_init.kernel_parameters.end(), kernel_parameters.begin());

          //append efficiencies to device vector
          thrust::copy(tmp_kernel_eff.begin(), tmp_kernel_eff.end(), kernel_parameters.begin()+n_user_params);

          // init kernel
          k_hall_davis_no_waals.resize(1, kernel_hall_davis_no_waals<real_t, n_t> (kernel_parameters.data(), detail::hall_davis_no_waals_r_max<real_t>()));
          p_kernel = (&(k_hall_davis_no_waals[0])).get();
          break;

        default:
          throw std::runtime_error("please supply a type of collision kernel to use");
      }
    }
  }
}
