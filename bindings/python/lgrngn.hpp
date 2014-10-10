// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

#include "error.hpp" 

#include <libcloudph++/lgrngn/factory.hpp>

namespace lgr = libcloudphxx::lgrngn;
namespace cmn = libcloudphxx::common;

namespace lgrngn
{
   lgr::particles_proto_t<real_t> *factory(
    const lgr::backend_t &backend,
    const lgr::opts_init_t<real_t> &opts_init
  ) {
    return lgr::factory(backend, opts_init);
  }

  bp::object outbuf(lgr::particles_proto_t<real_t> *arg)
  {
    int len = 1; // TODO: parcel assumption hardcoded
    return bp::object(bp::handle<>(PyBuffer_FromMemory(arg->outbuf(), len * sizeof(real_t)))); // TODO: this assumes Python 2 -> make it compatible with P3 or require P2 in CMake
  }

  void init(
    lgr::particles_proto_t<real_t> *arg,
    const bp::numeric::array &th,
    const bp::numeric::array &rv,
    const bp::numeric::array &rhod
  )
  {
    arg->init(
      np2ai(th),
      np2ai(rv),
      np2ai(rhod)
    );
    // TODO: 1D, 2D and 3D versions
  }

  void step_sync(
    lgr::particles_proto_t<real_t> *arg,
    const lgr::opts_t<real_t> &opts,
    const bp::numeric::array &th,
    const bp::numeric::array &rv,
    const bp::numeric::array &rhod
  )
  {
    lgr::arrinfo_t<real_t>
      np2ai_th(np2ai(th)),
      np2ai_rv(np2ai(rv));
    arg->step_sync(
      opts, 
      np2ai_th,
      np2ai_rv,
      np2ai(rhod)
    );
    // TODO: 1D, 2D and 3D versions
  }

  void set_dd(
    lgr::opts_init_t<real_t> *arg,
    const bp::dict &kappa_func
  )
  {
    for (int i = 0; i < len(kappa_func.keys()); ++i)
    {
      struct pyunary : cmn::unary_function<real_t> 
      {
        bp::object fun;

        pyunary(const bp::object &fun) : fun(fun) {}

	real_t funval(const real_t x) const
	{
	  return bp::extract<real_t>(fun(x)); 
	}
	
        pyunary *do_clone() const
        { 
          return new pyunary(*this); 
        }
      };

      boost::assign::ptr_map_insert<pyunary>(arg->dry_distros)(
        bp::extract<real_t>(kappa_func.keys()[i]), 
        pyunary(kappa_func.values()[i])
      );
    }
  }

  void get_dd(
    lgr::opts_init_t<real_t> *arg
  )
  {
    throw std::runtime_error("dry_distros does not feature a getter yet - TODO");
  }

  void set_cg(
    lgr::opts_t<real_t> *arg,
    const bp::dict &vec
  )
  {
    for (int i = 0; i < len(vec.keys()); ++i)
      arg->chem_gas.at(i) = bp::extract<real_t>(vec[i]);
  }

  bp::dict get_cg(
    lgr::opts_t<real_t> *arg
  )
  {
    bp::dict tmp;
    for (int i = 0; i < arg->chem_gas.size(); ++i)
      tmp[i] = arg->chem_gas.at(i);
    return tmp;
  }
};
