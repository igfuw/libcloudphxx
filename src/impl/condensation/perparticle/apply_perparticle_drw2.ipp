namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::apply_perparticle_drw2()
    {
      thrust_device::vector<real_t> &drw2 = drw2_gp->get(); 

      thrust::transform(
        rw2.begin(), rw2.end(),
        drw2.begin(),
        rw2.begin(),
        thrust::plus<real_t>()
      );
    }
  };
};