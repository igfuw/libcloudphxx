
namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sstp_percell_step(
    const int &step
    )
    {   
    if (sstp_cond == 1) return;

    namespace arg = thrust::placeholders;

    const int n = 3;
    thrust_device::vector<real_t>
        *scl[n] = { &rv,          &th,          &rhod        },
        *tmp[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh };

    for (int ix = 0; ix < (var_rho ? n : n-1); ++ix)
    {
        const real_t sstp = sstp_cond;
        if (step == 0)
        {
        thrust::transform(
            scl[ix]->begin(), scl[ix]->end(),
            tmp[ix]->begin(),
            tmp[ix]->begin(),
            arg::_1 - arg::_2
        );
        thrust::transform(
            scl[ix]->begin(), scl[ix]->end(),
            tmp[ix]->begin(),
            scl[ix]->begin(),
            arg::_1 - (sstp - 1) * arg::_2 / sstp
        );
        }
        else
        {
        thrust::transform(
            scl[ix]->begin(), scl[ix]->end(),
            tmp[ix]->begin(),
            scl[ix]->begin(),
            arg::_1 + arg::_2 / sstp
        );
        }
    }
    }
  };
};