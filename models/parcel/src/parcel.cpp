//
// Created by agnieszka on 08.09.23.
//

#include <libcloudph++/lgrngn/factory.hpp>
#include <iostream>

using real_t = float;

struct lognormal : libcloudphxx::common::unary_function<real_t> {
    const real_t mean_r, stdev, n_tot;

    explicit lognormal(
            const real_t mean_r, const real_t stdev, const real_t n_tot
    ) : mean_r(mean_r), stdev(stdev), n_tot(n_tot)
    {}

    real_t funval(const real_t lnr) const {
        return n_tot
        * real_t(exp(-pow((lnr - log(mean_r)), 2) / real_t(2) / pow(log(stdev),2)))
        / real_t(log(stdev))
        / real_t(sqrt(2*pi<real_t>()));
    }
};

struct bimodal : libcloudphxx::common::unary_function<real_t> {
    const lognormal lognormal1, lognormal2;

    explicit bimodal(
            const lognormal &lognormal1,
            const lognormal &lognormal2
    ) : lognormal1(lognormal1), lognormal2(lognormal2)
    {}

    real_t funval(const real_t lnr) const {
        return lognormal1(lnr) + lognormal2(lnr);
    }
};

int main(int argc, char** argv)
{
    std::unique_ptr<
        libcloudphxx::lgrngn::particles_proto_t<real_t>
    > prtcls;

    libcloudphxx::lgrngn::opts_init_t<real_t> params;

    real_t T = 283; // [K]
    real_t RH = .98;
    real_t p0 = 90000; // [Pa]

    real_t dv = 1; // [m3]
    real_t
      thd0 = 300, // [K]
      rv0 = 1e-3, // [kg/kg]
      rhod0 = 1; // [kg/m3]

    params.nx = params.ny = params.nz = 1;
    params.dx = params.dy = params.dz = pow(dv, real_t(1./3));
    params.dt = 1; // [s]
    params.sd_conc = 64; // n_sd
    params.n_sd_max = params.sd_conc;
    params.dry_distros.emplace(
        .666, // kappa
        std::make_shared<bimodal>(
                lognormal(1e-6, 1.2, 125e6),
                lognormal(1e-6, 1.2, 125e6)
        ) // n(ln(rd)) @ STP
    );
    params.coal_switch = false;
    params.sedi_switch = false;

    prtcls.reset(libcloudphxx::lgrngn::factory<real_t>(
        libcloudphxx::lgrngn::serial,
        params
    ));

    ptrdiff_t strides[] = {0, 0, 1};
    prtcls->init(
        libcloudphxx::lgrngn::arrinfo_t<real_t>(&thd0, strides),
        libcloudphxx::lgrngn::arrinfo_t<real_t>(&rv0, strides),
        libcloudphxx::lgrngn::arrinfo_t<real_t>(&rhod0, strides)
    );

    //prtcls->diag_all();
    prtcls->diag_dry_rng(0, 1e-6);
    prtcls->diag_sd_conc();
    std::cout << prtcls->outbuf()[0] << std::endl;
}