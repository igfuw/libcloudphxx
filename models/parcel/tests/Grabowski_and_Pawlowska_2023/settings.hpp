#pragma once

#include <cassert>
#include "spectral.hpp"

template <typename real_t>
class settings_t {
    const real_t
        T0 = 283, // [K]
        RH0 = .97,
        p0 = 90000; // [Pa]
public:
    const real_t
        z_max = 1000, // [m]
        vertical_velocity,
        dt = .01, // [s];
        kappa = .666;
    const int n_sd = 64;

    std::shared_ptr<bimodal<real_t>> n_ln_rd_stp;

    settings_t(
        const real_t vertical_velocity,
        const std::string aerosol = "pristine"
    ) : vertical_velocity(vertical_velocity) {
        if (aerosol == "pristine") {
            n_ln_rd_stp = std::make_shared<bimodal<real_t>>(
                    lognormal<real_t>(11*1e-9, 1.2, 125e6), //pristine
                    lognormal<real_t>(60*1e-9, 1.7, 65e6)
                    //lognormal(29*1e-9, 1.36, 160e6), //polluted
                    //lognormal(71*1e-9, 1.57, 380e6)
            ); // n(ln(rd)) @ STP
        }
        else assert(false);
    }

    auto n_steps() {
        auto n_steps = z_max / (vertical_velocity * dt);
        assert(n_steps == int(n_steps));
        return int(n_steps);
    }

    // TODO: compute from T0, p0, RH0
    real_t thd0() { return 292.; } // [K]
    real_t rv0() { return 0.008; } // [kg/kg]
    real_t rhod0() { return 1.1; } // [kg/m3]
};
