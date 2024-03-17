//
// Created by agnieszka on 19.09.23.
//

#pragma once
#include <netcdf>
#include <valarray>
#include <boost/range/combine.hpp>
#include "libcloud_hacks.hpp"
#include "settings.hpp"

template <typename real_t>
auto output_init(
        const int &n_sd,
        const int &n_t,
        const settings_t<real_t> &settings,
        const std::string &filename
) {
    auto nc = std::make_shared<netCDF::NcFile>(filename, netCDF::NcFile::replace);

    nc->putAtt("dt (s)", netCDF::ncFloat, settings.dt);
    nc->putAtt("vertical velocity (m s-1)", netCDF::ncFloat, settings.vertical_velocity);
    nc->putAtt("aerosol", settings.aerosol);
    nc->putAtt("init", settings.init);

    nc->addDim("step",n_t);
    nc->addDim("droplet_id",n_sd);

    nc->addVar("RH", "double", "step").putAtt("unit", "1");
    nc->addVar("T", "double", "step").putAtt("unit", "K");

    nc->addVar("wet radius squared", "double", std::vector<std::string>{"step", "droplet_id"}).putAtt("unit", "m^2");
    nc->addVar("dry radius cubed", "double", std::vector<std::string>{"step", "droplet_id"}).putAtt("unit", "m^3");
    nc->addVar("critical radius cubed", "double", std::vector<std::string>{"step", "droplet_id"}).putAtt("unit", "m^3");
    nc->addVar("kappa", "double", std::vector<std::string>{"step", "droplet_id"}).putAtt("unit", "1");

    return nc;
}

template <typename real_t>
auto save_scalar(
        const int i,
        libcloudphxx::lgrngn::particles_proto_t<real_t> &prtcls,
        const netCDF::NcFile &nc,
        const std::string &name
) {
    auto value = prtcls.outbuf()[0];
    nc.getVar(name).putVar(std::vector{size_t(i)}, value);
    return value;
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

template <libcloudphxx::lgrngn::backend_t backend, typename real_t>
const auto rw3_cr(
        libcloudphxx::lgrngn::particles_proto_t<real_t> &prtcls,
        const size_t &n_sd,
        const real_t &T
) {
    std::valarray<real_t> values(n_sd);
    for (size_t j = 0; const auto& [rd3, kpa]: boost::combine(
            impl<backend, real_t>(prtcls)->rd3,
            impl<backend, real_t>(prtcls)->kpa
    )) {
        values[j++] = libcloudphxx::common::kappa_koehler::rw3_cr(
                real_t(rd3) * si::cubic_meter,
                real_t(kpa) * si::dimensionless(),
                T * si::kelvin
        ) / si::cubic_meter;
    }
    return values;
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
    auto T = save_scalar(i, prtcls, nc, "T");

    save_vector(i, impl<backend, real_t>(prtcls)->rw2, nc, "wet radius squared");
    save_vector(i, impl<backend, real_t>(prtcls)->rd3, nc, "dry radius cubed");
    save_vector(i, impl<backend, real_t>(prtcls)->kpa, nc, "kappa");

    auto n_sd = nc.getDim("droplet_id").getSize();
    save_vector(i, rw3_cr<backend>(prtcls, n_sd, T), nc, "critical radius cubed");
}