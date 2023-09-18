//
// Created by agnieszka on 18.09.23.
//

#pragma once

#include "../../../../src/lib.hpp"
#include <thrust/system/cpp/vector.h>
namespace thrust_device = ::thrust::cpp;
#include "../../../../src/particles.tpp"

template <
    libcloudphxx::lgrngn::backend_t backend,
    typename real_t
>
auto impl(
    libcloudphxx::lgrngn::particles_proto_t<real_t> &prtcls
) {
    auto prtcls_with_backend = dynamic_cast<
            libcloudphxx::lgrngn::particles_t<real_t, backend>*
            >(&prtcls);
    return dynamic_cast<
            libcloudphxx::lgrngn::particles_t<real_t, backend>::impl*
            >(prtcls_with_backend->pimpl.get());
}
