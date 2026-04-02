
namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::acquire_arrays_for_perparticle_sstp(const bool cond, const bool depo)
    {        
      reset_guardp(sstp_dlt_rv_gp, tmp_device_real_part); 
      reset_guardp(sstp_dlt_th_gp, tmp_device_real_part); 
      reset_guardp(sstp_dlt_rhod_gp, tmp_device_real_part); 
      // if(opts_init.const_p)
          reset_guardp(sstp_dlt_p_gp, tmp_device_real_part); 

      if(!sstp_cond_exact_nomix_adaptive)
      {
        reset_guardp(Tp_gp, tmp_device_real_part);
        if(cond)
        {
          reset_guardp(rwX_gp, tmp_device_real_part);
          reset_guardp(drwX_gp, tmp_device_real_part);
        }
        if(depo)
        {
          reset_guardp(ice_mass_gp, tmp_device_real_part);
          reset_guardp(d_ice_mass_gp, tmp_device_real_part);
        }
      }

      if(opts_init.adaptive_sstp_cond)
      {
        reset_guardp(perparticle_sstp_cond_gp, tmp_device_n_part);
      }
    }
  };
};
