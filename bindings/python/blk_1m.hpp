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

namespace libcloudphxx
{
  namespace python
  {
    namespace b1m = libcloudphxx::blk_1m;

    namespace blk_1m
    {
      template <class arr_t>
      void adj_cellwise(
	const b1m::opts_t<typename arr_t::T_numtype>& opts,
	const bp_array &rhod,
	bp_array &th,
	bp_array &rv,
	bp_array &rc,
	bp_array &rr,
	const typename arr_t::T_numtype &dt
      )
      {
	arr_t 
	  np2bz_th(np2bz<arr_t>(th)), 
	  np2bz_rv(np2bz<arr_t>(rv)),
	  np2bz_rc(np2bz<arr_t>(rc)), 
	  np2bz_rr(np2bz<arr_t>(rr));
	b1m::adj_cellwise(
	  opts, 
	  np2bz<arr_t>(rhod), // since it is const, it may be a temporary object
	  np2bz_th, 
	  np2bz_rv, 
	  np2bz_rc, 
	  np2bz_rr, 
	  dt
	);
      }

      template <class arr_t>
      void adj_cellwise_constp(
	const b1m::opts_t<typename arr_t::T_numtype>& opts,
	const bp_array &rhod,
	const bp_array &p,
	bp_array &th,
	bp_array &rv,
	bp_array &rc,
	bp_array &rr,
	const typename arr_t::T_numtype &dt
      )
      {
	arr_t 
	  np2bz_th(np2bz<arr_t>(th)), 
	  np2bz_rv(np2bz<arr_t>(rv)),
	  np2bz_rc(np2bz<arr_t>(rc)), 
	  np2bz_rr(np2bz<arr_t>(rr));
	b1m::adj_cellwise_constp(
	  opts, 
	  np2bz<arr_t>(rhod), // since it is const, it may be a temporary object
	  np2bz<arr_t>(p), 
	  np2bz_th, 
	  np2bz_rv, 
	  np2bz_rc, 
	  np2bz_rr, 
	  dt
	);
      }

      template <class arr_t>
      void rhs_cellwise(
	const b1m::opts_t<typename arr_t::T_numtype> &opts,
	bp_array &dot_rc,
	bp_array &dot_rr,
	const bp_array &rc,
	const bp_array &rr
      ) 
      {
	arr_t
	  np2bz_dot_rc(np2bz<arr_t>(dot_rc)), 
	  np2bz_dot_rr(np2bz<arr_t>(dot_rr));
	b1m::rhs_cellwise(
	  opts,
	  np2bz_dot_rc,
	  np2bz_dot_rr,
	  np2bz<arr_t>(rc),
	  np2bz<arr_t>(rr)
	);
      } 

      template <class arr_t>
      typename arr_t::T_numtype rhs_columnwise(
	const b1m::opts_t<typename arr_t::T_numtype> &opts,
	bp_array &dot_rr,
	const bp_array &rhod,
	const bp_array &rr,
	const typename arr_t::T_numtype &dz
      ) {
	arr_t
	  np2bz_dot_rr(np2bz<arr_t>(dot_rr));
	return b1m::rhs_columnwise(
	  opts,
	  np2bz_dot_rr,
	  np2bz<arr_t>(rhod),
	  np2bz<arr_t>(rr),
	  dz
	);
      } 
    };
  };
};
