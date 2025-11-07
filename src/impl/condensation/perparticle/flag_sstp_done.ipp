namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::flag_sstp_done(const int step_plus_one)
    {
      namespace arg = thrust::placeholders;

      auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get(); //note the use of unconverged_mask!
      auto &sstp_cond = perparticle_sstp_cond_gp->get();

      // debug::print(sstp_cond);

      thrust::transform_if(        
        sstp_cond.begin(),
        sstp_cond.end(),
        unconverged_mask.begin(),
        unconverged_mask.begin(),
        arg::_1 != step_plus_one,
        cuda::std::identity()
      );

      // debug::print(unconverged_mask);

      // std::cerr << "particles done in step " << step_plus_one << ": " 
      //           << thrust::count(unconverged_mask.begin(), unconverged_mask.end(), false) 
      //           << "/" << n_part << std::endl;
    }
  };
};