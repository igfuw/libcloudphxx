// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

#include "error.hpp"

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

namespace b2m = libcloudphxx::blk_2m;

namespace blk_2m
{
  void rhs_cellwise(
    const b2m::opts_t<real_t> &opts,
    bp::numeric::array &dot_th,
    bp::numeric::array &dot_rv,
    bp::numeric::array &dot_rc,
    bp::numeric::array &dot_nc,
    bp::numeric::array &dot_rr,
    bp::numeric::array &dot_nr,
    const bp::numeric::array &rhod,
    const bp::numeric::array &th,
    const bp::numeric::array &rv,
    const bp::numeric::array &rc,
    const bp::numeric::array &nc,
    const bp::numeric::array &rr,
    const bp::numeric::array &nr,
    const real_t &dt
  ) 
  {
    arr_t
      np2bz_dot_th(np2bz(dot_th)), 
      np2bz_dot_rv(np2bz(dot_rv)), 
      np2bz_dot_rc(np2bz(dot_rc)), 
      np2bz_dot_nc(np2bz(dot_nc)), 
      np2bz_dot_rr(np2bz(dot_rr)), 
      np2bz_dot_nr(np2bz(dot_nr));
    b2m::rhs_cellwise(
      opts,
      np2bz_dot_th,
      np2bz_dot_rv,
      np2bz_dot_rc,
      np2bz_dot_nc,
      np2bz_dot_rr,
      np2bz_dot_nr,
      np2bz(rhod),
      np2bz(th),
      np2bz(rv),
      np2bz(rc),
      np2bz(nc),
      np2bz(rr),
      np2bz(nr),
      dt
    );
  } 

  real_t rhs_columnwise(
    const b2m::opts_t<real_t> &opts,
    bp::numeric::array &dot_rr,
    bp::numeric::array &dot_nr,
    const bp::numeric::array &rhod,
    const bp::numeric::array &rr,
    const bp::numeric::array &nr,
    const real_t &dt,
    const real_t &dz
  ) {
    arr_t
      np2bz_dot_rr(np2bz(dot_rr)), 
      np2bz_dot_nr(np2bz(dot_nr));
    return b2m::rhs_columnwise(
      opts,
      np2bz_dot_rr,
      np2bz_dot_nr,
      np2bz(rhod),
      np2bz(rr),
      np2bz(nr),
      dt,
      dz
    );
  } 

  void set_dd(
    b2m::opts_t<real_t> *arg,
    const bp::list &modes
  )
  {
    for (int i = 0; i < len(modes); ++i)
    {
      arg->dry_distros.push_back({
        .mean_rd = bp::extract<real_t>(modes[i]["mean_rd"]),
        .sdev_rd = bp::extract<real_t>(modes[i]["sdev_rd"]),
        .N_stp   = bp::extract<real_t>(modes[i]["N_stp"]),
        .chem_b  = bp::extract<real_t>(modes[i]["chem_b"])
      });
    }
  }

  void get_dd(
    b2m::opts_t<real_t> *arg
  )
  {
    throw std::runtime_error("dry_distros does not feature a getter yet - TODO");
  }
};
