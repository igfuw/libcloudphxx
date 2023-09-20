//
// Created by agnieszka on 19.09.23.
//

#pragma once

#include <libcloudph++/common/earth.hpp>


template <typename real_t>
void step_hydrostatic(
    const real_t &dz_arg, const real_t &thd_arg, const real_t &rv_arg,
    real_t &rhod_arg
) {
    auto dz = dz_arg * si::meter;
    auto thd = thd_arg * si::kelvin;
    auto rv = rv_arg * si::dimensionless();
    auto rhod = rhod_arg * si::kilogram_per_cubic_meter;
    auto p0 = libcloudphxx::common::theta_std::p_1000<real_t>();
    auto
            R = libcloudphxx::common::moist_air::R(rv),
            cp = libcloudphxx::common::moist_air::c_p(rv);
    auto K = R/cp;
    auto g = libcloudphxx::common::earth::g<real_t>();

    //using boost::units::pow;
    auto rho = rhod/(1-rv);
    auto T = libcloudphxx::common::theta_dry::T(thd, rhod);
    auto p = libcloudphxx::common::theta_dry::p(rhod, rv, T);
//    auto d_rho_dz = pow(rho,(2*K-1)/(K-1))*(K-1)*g*pow(R*thd/real_t(pow(p0,K)),1/(K-1));
    auto d_rho_dz = g / T * rho * (R / cp - 1) / R;
    rho += d_rho_dz * dz;
    rhod=rho*(1-rv);

    rhod_arg = rhod / si::kilogram_per_cubic_meter;
}
