namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::apply_perparticle_cond_change_to_percell_rv_and_th()
    {
      if(opts_init.sstp_cond_mix)
      {
        update_state(rv, sstp_tmp_rv);
        update_state(th, sstp_tmp_th);
      }
      else
      {
        rw_mom3_post_change();  
        update_th_rv();
      }
    }
  };
};