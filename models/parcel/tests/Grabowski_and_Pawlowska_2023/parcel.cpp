//
// Created by agnieszka on 08.09.23.
//

#include "libcloud_hacks.hpp"

#include <libcloudph++/lgrngn/factory.hpp>
#include <iostream>
#include <ranges>
#include "settings.hpp"


#include <gnuplot-iostream.h>



using real_t = float;
const auto backend = libcloudphxx::lgrngn::serial;



void output(libcloudphxx::lgrngn::particles_proto_t<real_t> &prtcls) {
    //prtcls->diag_all();
    prtcls.diag_dry_rng(0, 1e-6);
    prtcls.diag_sd_conc();
    auto n_sd = prtcls.outbuf()[0];

    prtcls.diag_RH();
    auto RH = prtcls.outbuf()[0];


    Gnuplot gp;
    gp << "set term svg size " << 500 << "," << 500 << "\n";
    gp << "set output 'test.svg'"<<"\n";
    gp << "plot '-'" << "\n";
    gp.send1d(impl<backend, real_t>(prtcls)->rd3);

    for (const auto &rd3 : impl<backend, real_t>(prtcls)->rw2) // rd3, n, kpa
        std::cout << rd3 << std::endl;

    //std::cout << n_sd << "\t" << RH << std::endl;


}

template <typename real_t>
auto arrinfo(real_t &var) {
    static const ptrdiff_t strides[] = {0, 0, 1};
    return libcloudphxx::lgrngn::arrinfo_t(&var, strides);
}

template <typename real_t>
void step_hydrostatic(const real_t &dz, const real_t &thd, const real_t &rv, real_t &rhod) {
    real_t Rv= 461.5, Rd=287.04, cpv=1875, cpd=1004.67, g=9.81, p0=100000;
    real_t R=rv*Rv+(1-rv)*Rd, cp=rv*cpv+(1-rv)*cpd, K=R/cp;

    real_t rho = rhod/(1-rv);
    real_t d_rho_dz = pow(rho,(2*K-1)/(K-1))*(K-1)*g*pow(R*thd/pow(p0,K),1/(K-1));
    rho += d_rho_dz * dz;
    rhod=rho*(1-rv);
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
        backend,
        params
    ));



    real_t
        thd = settings.thd0(),
        rv = settings.rv0(),
        rhod = settings.rhod0();

    prtcls->init(arrinfo(thd), arrinfo(rv), arrinfo(rhod));
    for (auto i : std::views::iota(0, settings.n_steps() + 1)) {
        if (i != 0) {
            step_hydrostatic(settings.dt*settings.vertical_velocity, thd, rv, rhod);
            prtcls->step_sync(opts, arrinfo(thd), arrinfo(rv), arrinfo(rhod));
            prtcls->step_async(opts);
        }
        output(*prtcls);
    }
}
