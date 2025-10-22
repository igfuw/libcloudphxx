
namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::acquire_arrays_for_perparticle_sstp()
    {        
    reset_guardp(sstp_dlt_rv_gp, tmp_device_real_part); 
    reset_guardp(sstp_dlt_th_gp, tmp_device_real_part); 
    reset_guardp(sstp_dlt_rhod_gp, tmp_device_real_part); 
    if(opts_init.const_p)
        reset_guardp(sstp_dlt_p_gp, tmp_device_real_part); 

    reset_guardp(rw3_gp, tmp_device_real_part);
    reset_guardp(drw3_gp, tmp_device_real_part);
    reset_guardp(Tp_gp, tmp_device_real_part);

    if(opts_init.adaptive_sstp_cond)
        reset_guardp(sstp_cond_perparticle_gp, tmp_device_n_part);
    }
  };
};