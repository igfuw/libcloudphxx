namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::store_unconverged_perparticle_drw2_as_old(
      const thrust_device::vector<real_t> &drw2,
      thrust_device::vector<real_t> &drw2_old
    ) { 
      const auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();

      thrust::transform_if(
        drw2.begin(), drw2.end(),
        unconverged_mask.begin(),
        drw2_old.begin(),
        cuda::std::identity(),
        cuda::std::identity()
      );
    }
  };
};