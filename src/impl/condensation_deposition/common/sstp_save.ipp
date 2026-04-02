
namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sstp_save()
    {
    if (!allow_sstp_cond) return;

    const int n = 4;
    thrust_device::vector<real_t>
        *fr[n] = { &rv,          &th,          &rhod,        &p          },
        *to[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh, &sstp_tmp_p };
    for (int ix = 0; ix < ( (opts_init.const_p && opts_init.exact_sstp_cond) ? n : n-1); ++ix) // TODO: var_rho
    {
        if(opts_init.exact_sstp_cond) // per-particle version
        thrust::copy(
            thrust::make_permutation_iterator(fr[ix]->begin(), ijk.begin()),
            thrust::make_permutation_iterator(fr[ix]->begin(), ijk.end()),
            to[ix]->begin()
        ); 
        else // per-cell version
        thrust::copy(
            fr[ix]->begin(),
            fr[ix]->end(),
            to[ix]->begin()
        ); 
    }
    }
  };
};