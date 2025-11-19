
namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::release_arrays_for_perparticle_sstp()
    {        
      sstp_dlt_rv_gp.reset();
      sstp_dlt_th_gp.reset();
      sstp_dlt_rhod_gp.reset();
      // if(opts_init.const_p)
          sstp_dlt_p_gp.reset();

      rwX_gp.reset();
      drwX_gp.reset();
      // drw2_gp.reset();
      // drw3_gp.reset();
      Tp_gp.reset();

      if(opts_init.adaptive_sstp_cond)
      {
          perparticle_sstp_cond_gp.reset();
      }
    }
  };
};