#pragma once
#include <libcloudph++/common/unary_function.hpp>

template <typename real_t>
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

template <typename real_t>
struct bimodal : libcloudphxx::common::unary_function<real_t> {
    const lognormal<real_t> lognormal1, lognormal2;

    explicit bimodal(
            const lognormal<real_t> &lognormal1,
            const lognormal<real_t> &lognormal2
    ) : lognormal1(lognormal1), lognormal2(lognormal2)
    {}

    real_t funval(const real_t lnr) const {
        return lognormal1(lnr) + lognormal2(lnr);
    }
};