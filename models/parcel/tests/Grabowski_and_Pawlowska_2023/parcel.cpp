//
// Created by agnieszka on 08.09.23.
//

#include <ranges>

#include "libcloud_hacks.hpp"
#include "settings.hpp"
#include "hydrostatics.hpp"
#include "output.hpp"
#include "command_line_options.hpp"

#include <libcloudph++/lgrngn/factory.hpp>

using real_t = float;
const auto backend = libcloudphxx::lgrngn::serial;


int main(int arg_count, char** arg_values)
{
    auto vm = parse_options(arg_count, arg_values);
    settings_t<real_t> settings(
        vm["vertical_velocity"].as<real_t>(),
        vm["aerosol"].as<std::string>(),
        vm["dt"].as<real_t>()
    );

    libcloudphxx::lgrngn::opts_init_t<real_t> params;
    libcloudphxx::lgrngn::opts_t<real_t> opts;

    real_t dv = 1; // [m3]
    params.nx = params.ny = params.nz = 1;
    params.dx = params.dy = params.dz = pow(dv, real_t(1./3));
    params.dt = settings.dt;
    params.sd_conc = params.n_sd_max = settings.n_sd;
    params.dry_distros.emplace(settings.kappa, settings.n_ln_rd_stp);
    params.sstp_cond = 10;

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

    auto filename=vm["filename"].as<std::string>();
    auto range_i = std::views::iota(0, settings.n_steps() + 1);
    auto nc = output_init(settings.n_sd, range_i.size(), settings, filename);

    prtcls->init(arrinfo(thd), arrinfo(rv), arrinfo(rhod));
    {
        for (auto i: range_i) {
            if (i != 0) {
                step_hydrostatic(settings.dt * settings.vertical_velocity, thd, rv, rhod);
                prtcls->step_sync(opts, arrinfo(thd), arrinfo(rv), arrinfo(rhod));
                prtcls->step_async(opts);
            }
            output_step<backend>(i, *prtcls, *nc);
        }
    }
}
