//
// Created by agnieszka on 08.09.23.
//

#include <libcloudph++/lgrngn/factory.hpp>
#include <iostream>
#include <ranges>
#include "settings.hpp"

using real_t = float;


void output(libcloudphxx::lgrngn::particles_proto_t<real_t> &prtcls) {
    //prtcls->diag_all();
    prtcls.diag_dry_rng(0, 1e-6);
    prtcls.diag_sd_conc();
    auto n_sd = prtcls.outbuf()[0];

    prtcls.diag_RH();
    auto RH = prtcls.outbuf()[0];

    std::cout << n_sd << "\t" << RH << std::endl;
}

template <typename real_t>
auto ai(real_t &var) {
    static const ptrdiff_t strides[] = {0, 0, 1};
    return libcloudphxx::lgrngn::arrinfo_t(&var, strides);
}

template <typename real_t>
void step_hydrostatic(const real_t &dt, const real_t &thd, const real_t &rv, real_t &rhod) {
    real_t d_rhod_dt = -.0001; // kg/m3/s
    rhod += d_rhod_dt * dt;
}

int main(int argc, char** argv)
{
    settings_t<real_t> settings(1.);

    libcloudphxx::lgrngn::opts_init_t<real_t> params;
    libcloudphxx::lgrngn::opts_t<real_t> opts;

    real_t dv = 1; // [m3]
    params.nx = params.ny = params.nz = 1;
    params.dx = params.dy = params.dz = pow(dv, real_t(1./3));
    params.dt = settings.dt;
    params.sd_conc = params.n_sd_max = settings.n_sd;
    params.dry_distros.emplace(settings.kappa, settings.n_ln_rd_stp);

    params.coal_switch = opts.coal = false;
    params.sedi_switch = opts.sedi = false;

    std::unique_ptr<
        libcloudphxx::lgrngn::particles_proto_t<real_t>
    > prtcls(libcloudphxx::lgrngn::factory<real_t>(
        libcloudphxx::lgrngn::serial,
        params
    ));

    real_t
        thd = settings.thd0(),
        rv = settings.rv0(),
        rhod = settings.rhod0();

    prtcls->init(ai(thd), ai(rv), ai(rhod));
    for (auto i : std::views::iota(0, settings.n_steps() + 1)) {
        if (i != 0) {
            step_hydrostatic(settings.dt, thd, rv, rhod);
            prtcls->step_sync(opts, ai(thd), ai(rv), ai(rhod));
            prtcls->step_async(opts);
        }
        output(*prtcls);
    }
}
