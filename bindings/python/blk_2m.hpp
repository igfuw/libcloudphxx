// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

#include "error.hpp"

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

namespace libcloudphxx
{
  namespace python
  {
    namespace b2m = libcloudphxx::blk_2m;

    namespace blk_2m
    {
      template <typename arr_t>
      void rhs_cellwise(
	const b2m::opts_t<typename arr_t::T_numtype> &opts,
	bp_array &dot_th,
	bp_array &dot_rv,
	bp_array &dot_rc,
	bp_array &dot_nc,
	bp_array &dot_rr,
	bp_array &dot_nr,
	const bp_array &rhod,
	const bp_array &th,
	const bp_array &rv,
	const bp_array &rc,
	const bp_array &nc,
	const bp_array &rr,
	const bp_array &nr,
	const typename arr_t::T_numtype &dt
      )
      {
	arr_t
	  np2bz_dot_th(np2bz<arr_t>(dot_th)),
	  np2bz_dot_rv(np2bz<arr_t>(dot_rv)),
	  np2bz_dot_rc(np2bz<arr_t>(dot_rc)),
	  np2bz_dot_nc(np2bz<arr_t>(dot_nc)),
	  np2bz_dot_rr(np2bz<arr_t>(dot_rr)),
	  np2bz_dot_nr(np2bz<arr_t>(dot_nr)),
	  np2bz_nc(np2bz<arr_t>(nc)),
	  np2bz_rc(np2bz<arr_t>(rc)),
	  np2bz_nr(np2bz<arr_t>(nr)),
	  np2bz_rr(np2bz<arr_t>(rr));
	b2m::rhs_cellwise(
	  opts,
	  np2bz_dot_th,
	  np2bz_dot_rv,
	  np2bz_dot_rc,
	  np2bz_dot_nc,
	  np2bz_dot_rr,
	  np2bz_dot_nr,
	  np2bz<arr_t>(rhod),
	  np2bz<arr_t>(th),
	  np2bz<arr_t>(rv),
	  np2bz_rc,
	  np2bz_nc,
	  np2bz_rr,
	  np2bz_nr,
	  dt
	);
      }

      template <typename arr_t>
      typename arr_t::T_numtype rhs_columnwise(
	const b2m::opts_t<typename arr_t::T_numtype> &opts,
	bp_array &dot_rr,
	bp_array &dot_nr,
	const bp_array &rhod,
	const bp_array &rr,
	const bp_array &nr,
	const typename arr_t::T_numtype &dt,
	const typename arr_t::T_numtype &dz
      ) {
	arr_t
	  np2bz_dot_rr(np2bz<arr_t>(dot_rr)),
	  np2bz_dot_nr(np2bz<arr_t>(dot_nr));
	return b2m::rhs_columnwise(
	  opts,
	  np2bz_dot_rr,
	  np2bz_dot_nr,
	  np2bz<arr_t>(rhod),
	  np2bz<arr_t>(rr),
	  np2bz<arr_t>(nr),
	  dt,
	  dz
	);
      }

      template <typename real_t>
      void set_dd(
	b2m::opts_t<real_t> *arg,
	const bp::list &modes
      )
      {
	for (int i = 0; i < len(modes); ++i)
	{
	  arg->dry_distros.push_back({
	     bp::extract<real_t>(modes[i]["mean_rd"]),
	     bp::extract<real_t>(modes[i]["sdev_rd"]),
	     bp::extract<real_t>(modes[i]["N_stp"]),
	     bp::extract<real_t>(modes[i]["chem_b"])
	  });
	}
      }

      template <typename real_t>
      void get_dd(
	b2m::opts_t<real_t> *arg
      )
      {
	throw std::runtime_error("dry_distros does not feature a getter yet - TODO");
      }
    };
  };
};
