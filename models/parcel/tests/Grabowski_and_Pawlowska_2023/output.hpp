//
// Created by agnieszka on 19.09.23.
//

#pragma once
#include <netcdf>
#include "libcloud_hacks.hpp"
#include "settings.hpp"

template <typename real_t>
auto output_init(
    const int &n_sd,
    const int &n_t,
    const settings_t<real_t> &settings
) {
    auto nc = std::make_shared<netCDF::NcFile>("test.nc", netCDF::NcFile::replace);

    nc->putAtt("dt (s)", netCDF::ncFloat, settings.dt);
    nc->putAtt("vertical velocity (m s-1)", netCDF::ncFloat, settings.vertical_velocity);
    nc->putAtt("aerosol", settings.aerosol);

    nc->addDim("step",n_t);
    nc->addDim("droplet_id",n_sd);

    nc->addVar("RH", "float", "step").putAtt("unit", "1");
    nc->addVar("T", "float", "step").putAtt("unit", "K");

    nc->addVar("wet radius squared", "float", std::vector<std::string>{"step", "droplet_id"}).putAtt("unit", "m^2");
    nc->addVar("dry radius cubed", "float", std::vector<std::string>{"step", "droplet_id"}).putAtt("unit", "m^3");
    nc->addVar("kappa", "float", std::vector<std::string>{"step", "droplet_id"}).putAtt("unit", "1");

    return nc;
}


template <typename real_t>
void save_scalar(
    const int i,
    libcloudphxx::lgrngn::particles_proto_t<real_t> &prtcls,
    const netCDF::NcFile &nc,
    const std::string &name
) {
    nc.getVar(name).putVar(std::vector{size_t(i)}, prtcls.outbuf()[0]);
}


template <typename iter_t>
void save_vector(
        const int i,
        iter_t &iter,
        const netCDF::NcFile &nc,
        const std::string &name
)
{
    for (size_t j=0; const auto &item : iter)
        nc.getVar(name).putVar(std::vector{size_t(i), j++}, item);
}

template <
        libcloudphxx::lgrngn::backend_t backend,
        typename real_t
>
void output_step(
        const int i,
        libcloudphxx::lgrngn::particles_proto_t<real_t> &prtcls,
        const netCDF::NcFile &nc
) {
    prtcls.diag_RH();
    save_scalar(i, prtcls, nc, "RH");

    prtcls.diag_temperature();
    save_scalar(i, prtcls, nc, "T");

    save_vector(i, impl<backend, real_t>(prtcls)->rw2, nc, "wet radius squared");
    save_vector(i, impl<backend, real_t>(prtcls)->rd3, nc, "dry radius cubed");
    save_vector(i, impl<backend, real_t>(prtcls)->kpa, nc, "kappa");
}
