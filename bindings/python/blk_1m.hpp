// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#include "error.hpp"

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_columnwise.hpp>

namespace b1m = libcloudphxx::blk_1m;

namespace blk_1m
{
  void adj_cellwise(
    const b1m::opts_t<real_t>& opts,
    const bp::numeric::array &rhod,
    bp::numeric::array &th,
    bp::numeric::array &rv,
    bp::numeric::array &rc,
    bp::numeric::array &rr,
    const real_t &dt
  )
  {
    arr_t 
      np2bz_th(np2bz(th)), 
      np2bz_rv(np2bz(rv)),
      np2bz_rc(np2bz(rc)), 
      np2bz_rr(np2bz(rr));
    b1m::adj_cellwise(
      opts, 
      np2bz(rhod), // since it is const, it may be a temporary object
      np2bz_th, 
      np2bz_rv, 
      np2bz_rc, 
      np2bz_rr, 
      dt
    );
  }

  void rhs_cellwise(
    const b1m::opts_t<real_t> &opts,
    bp::numeric::array &dot_rc,
    bp::numeric::array &dot_rr,
    const bp::numeric::array &rc,
    const bp::numeric::array &rr
  ) 
  {
    arr_t
      np2bz_dot_rc(np2bz(dot_rc)), 
      np2bz_dot_rr(np2bz(dot_rr));
    b1m::rhs_cellwise(
      opts,
      np2bz_dot_rc,
      np2bz_dot_rr,
      np2bz(rc),
      np2bz(rr)
    );
  } 

  real_t rhs_columnwise(
    const b1m::opts_t<real_t> &opts,
    bp::numeric::array &dot_rr,
    const bp::numeric::array &rhod,
    const bp::numeric::array &rr,
    const real_t &dz
  ) {
    arr_t
      np2bz_dot_rr(np2bz(dot_rr));
    return b1m::rhs_columnwise(
      opts,
      np2bz_dot_rr,
      np2bz(rhod),
      np2bz(rr),
      dz
    );
  } 
};
