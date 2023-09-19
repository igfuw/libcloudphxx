//
// Created by agnieszka on 19.09.23.
//

#pragma once
#include <fstream>
#include "libcloud_hacks.hpp"

template <
        libcloudphxx::lgrngn::backend_t backend,
        typename real_t
>
void output(libcloudphxx::lgrngn::particles_proto_t<real_t> &prtcls, std::ofstream &myfile) {
    //prtcls->diag_all();
    prtcls.diag_dry_rng(0, 1e-6);
    prtcls.diag_sd_conc();
    auto n_sd = prtcls.outbuf()[0];

    prtcls.diag_RH();
    auto RH = prtcls.outbuf()[0];


    //for (const auto &r2 : impl<backend, real_t>(prtcls)->rw2) // rd3, n, kpa
    //   std::cout << r2 << std::endl;

    for (const auto &r2 : impl<backend, real_t>(prtcls)->rw2) { // rd3, n, kpa
        myfile << r2 << ",";
    }
    myfile<<std::endl;



    //std::cout << n_sd << "\t" << RH << std::endl;


}
