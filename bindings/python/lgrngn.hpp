// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

#include "error.hpp" 

#include <boost/assign/ptr_map_inserter.hpp>  // for 'ptr_map_insert()'

#include <libcloudph++/lgrngn/factory.hpp>

namespace libcloudphxx
{
  namespace python
  {
    namespace lgr = libcloudphxx::lgrngn;
    namespace cmn = libcloudphxx::common;

    namespace lgrngn
    {
      template <typename real_t>
      lgr::particles_proto_t<real_t> *factory(
	const lgr::backend_t &backend,
	const lgr::opts_init_t<real_t> &opts_init
      ) {
	return lgr::factory(backend, opts_init);
      }

      template <typename real_t>
      bp::object outbuf(
        lgr::particles_proto_t<real_t> *arg
      ) {
	return bp::object(bp::handle<>(PyBuffer_FromMemory(
          arg->outbuf(), 
          sizeof(real_t)
          * std::max(1, arg->opts_init->nx) 
          * std::max(1, arg->opts_init->ny) 
          * std::max(1, arg->opts_init->nz) 
        ))); // TODO: this assumes Python 2 -> make it compatible with P3 or require P2 in CMake
      }

      template <typename real_t>
      const std::array<int, 3> sz(
        const lgr::particles_proto_t<real_t> &arg
      ) {
        return std::array<int, 3>({
          arg.opts_init->nx,
          arg.opts_init->ny,
          arg.opts_init->nz
        });
      }

      // 0D or non-kinematic version
      template <typename real_t>
      void init_3arg(
	lgr::particles_proto_t<real_t> *arg,
	const bp::numeric::array &th,
	const bp::numeric::array &rv,
	const bp::numeric::array &rhod
      )
      {
	arg->init(
	  np2ai<real_t>(th,   sz(*arg)),
	  np2ai<real_t>(rv,   sz(*arg)),
	  np2ai<real_t>(rhod, sz(*arg))
	);
      }

      // 2D kinematic version
      template <typename real_t>
      void init_5arg(
	lgr::particles_proto_t<real_t> *arg,
	const bp::numeric::array &th,
	const bp::numeric::array &rv,
	const bp::numeric::array &rhod,
        const bp::numeric::array &rhod_Cx,
        const bp::numeric::array &rhod_Cz
      )
      {
	arg->init(
	  np2ai<real_t>(th,      sz(*arg)),
	  np2ai<real_t>(rv,      sz(*arg)),
	  np2ai<real_t>(rhod,    sz(*arg)),
          np2ai<real_t>(rhod_Cx, sz(*arg)),
          np2ai<real_t>(rhod_Cz, sz(*arg))
	);
      }

      // TODO: 1D & 3D kinematic versions

      template <typename real_t>
      void step_sync(
	lgr::particles_proto_t<real_t> *arg,
	const lgr::opts_t<real_t> &opts,
	const bp::numeric::array &th,
	const bp::numeric::array &rv,
	const bp::numeric::array &rhod
        // TODO: courant fields
      )
      {
	lgr::arrinfo_t<real_t>
	  np2ai_th(np2ai<real_t>(th, sz(*arg))),
	  np2ai_rv(np2ai<real_t>(rv, sz(*arg)));
	arg->step_sync(
	  opts, 
	  np2ai_th,
	  np2ai_rv,
	  np2ai<real_t>(rhod, sz(*arg))
	);
      }

      template <typename real_t>
      const lgr::opts_init_t<real_t> get_oi(
        lgr::particles_proto_t<real_t> *arg
      )
      {
        return *arg->opts_init;
      }

      template <typename real_t>
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

      template <typename real_t>
      void get_dd(
	lgr::opts_init_t<real_t> *arg
      )
      {
	throw std::runtime_error("dry_distros does not feature a getter yet - TODO");
      }

      template <typename real_t>
      void set_cg(
	lgr::opts_t<real_t> *arg,
	const bp::dict &vec
      )
      {
	for (int i = 0; i < len(vec.keys()); ++i)
	  arg->chem_gas.at(i) = bp::extract<real_t>(vec[i]);
      }

      template <typename real_t>
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
  };
};
