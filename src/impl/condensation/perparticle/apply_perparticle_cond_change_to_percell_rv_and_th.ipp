namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::apply_perparticle_cond_change_to_percell_rv_and_th()
    {
      if(opts_init.sstp_cond_mix)
      {
        // with mixing, all SDs in a cell should arrive at the same sstp_tmp_rv/th;
        // however, there may be some small divergences, as we see some differences in results
        // of the physics/lgrngn_cond_sstp.py test between per-cell and per-particle substepping with sstp_cond_mix=True;
        // these differences started to appear after changing the way how drw3 is applied to rv and th in per-particle substepping
        // (commit e78392573eedc6acd74861d5a5ae57a82d2edff0). This commit did not affect per-particle substepping without mixing.
        // Discussed divergences can be decreased (but not completly removed) by updating cell rv/th based on mom3, as done without mixing.
        // Nevertheless, these divergences are small enough to be acceptable in most applications.
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