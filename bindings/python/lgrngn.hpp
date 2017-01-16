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
      namespace detail
      {
        template<class real_t>
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
      };

      template <typename real_t>
      lgr::particles_proto_t<real_t> *factory(
	const lgr::backend_t &backend,
	lgr::opts_init_t<real_t> opts_init
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

      // 
      template <typename real_t>
      void init(
	lgr::particles_proto_t<real_t> *arg,
	const bp::numeric::array &th,
	const bp::numeric::array &rv,
	const bp::numeric::array &rhod,
        const bp::numeric::array &Cx,
        const bp::numeric::array &Cy,
        const bp::numeric::array &Cz,
        const bp::dict &ambient_chem
      )
      {
        typedef std::map<enum lgr::chem_species_t, const lgr::arrinfo_t<real_t> > map_t;
        map_t map;

        for (int i = 0; i < len(ambient_chem.keys()); ++i)
          map.insert(typename map_t::value_type(
            bp::extract<enum lgr::chem_species_t>(ambient_chem.keys()[i]),
            np2ai<real_t>(bp::extract<bp::numeric::array>(ambient_chem.values()[i]), sz(*arg))
          ));

	arg->init(
	  np2ai<real_t>(th,      sz(*arg)),
	  np2ai<real_t>(rv,      sz(*arg)),
	  np2ai<real_t>(rhod,    sz(*arg)),
          np2ai<real_t>(Cx,      sz(*arg)),
          np2ai<real_t>(Cy,      sz(*arg)),
          np2ai<real_t>(Cz,      sz(*arg)),
          map // ambient_chem
	);
      }

      // 
      template <typename real_t>
      void step_sync(
	lgr::particles_proto_t<real_t> *arg,
	const lgr::opts_t<real_t> &opts,
	const bp::numeric::array &th,
	const bp::numeric::array &rv,
	const bp::numeric::array &rhod,
	const bp::numeric::array &Cx,
	const bp::numeric::array &Cy,
	const bp::numeric::array &Cz,
        bp::dict &ambient_chem
      )
      {
        typedef std::map<enum lgr::chem_species_t, lgr::arrinfo_t<real_t> > map_t;
        map_t map;

        for (int i = 0; i < len(ambient_chem.keys()); ++i)
          map.insert(typename map_t::value_type(
            bp::extract<enum lgr::chem_species_t>(ambient_chem.keys()[i]),
            np2ai<real_t>(bp::extract<bp::numeric::array>(ambient_chem.values()[i]), sz(*arg))
          ));

	lgr::arrinfo_t<real_t>
	  np2ai_th(np2ai<real_t>(th, sz(*arg))),
	  np2ai_rv(np2ai<real_t>(rv, sz(*arg)));
	arg->step_sync(
	  opts, 
	  np2ai_th,
	  np2ai_rv,
	  np2ai<real_t>(rhod, sz(*arg)),
	  np2ai<real_t>(Cx, sz(*arg)),
	  np2ai<real_t>(Cy, sz(*arg)),
	  np2ai<real_t>(Cz, sz(*arg)),
          map
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
      void set_dd( // dry_distro
	lgr::opts_init_t<real_t> *arg,
	const bp::dict &kappa_func)
      {
        arg->dry_distros.clear();
        for (int i = 0; i < len(kappa_func.keys()); ++i)
          boost::assign::ptr_map_insert<detail::pyunary<real_t>>(arg->dry_distros)(
            bp::extract<real_t>(kappa_func.keys()[i]), 
            detail::pyunary<real_t>(kappa_func.values()[i])
          );
      }

      template <typename real_t>
      void set_sdd( // src_dry_distro
	lgr::opts_init_t<real_t> *arg,
	const bp::dict &kappa_func)
      {
        arg->src_dry_distros.clear();
        for (int i = 0; i < len(kappa_func.keys()); ++i)
          boost::assign::ptr_map_insert<detail::pyunary<real_t>>(arg->src_dry_distros)(
            bp::extract<real_t>(kappa_func.keys()[i]), 
            detail::pyunary<real_t>(kappa_func.values()[i])
          );
      }

      template <typename real_t>
      void set_ds( // dry_sizes
	lgr::opts_init_t<real_t> *arg,
	const bp::dict &kappa_func
      )
      {
        // TODO: loop over kappas (right now only one possible)
        const bp::dict size_conc = bp::extract<bp::dict>(kappa_func.values()[0]);
        std::map<real_t, real_t> size_conc_map;

        // turn the size-conc dict into a size-conc map
	for (int i = 0; i < len(size_conc.keys()); ++i)
        {
          size_conc_map[bp::extract<real_t>(size_conc.keys()[i])] = bp::extract<real_t>(size_conc.values()[i]);
          std::cout << "bindings map: " << bp::extract<real_t>(size_conc.keys()[i]) << " " <<
            size_conc_map[bp::extract<real_t>(size_conc.keys()[i])] << std::endl;;
        }
        const real_t kappa = bp::extract<real_t>(kappa_func.keys()[0]);
        arg->dry_sizes.clear();
        arg->dry_sizes[kappa] = size_conc_map;
      }

      template <typename real_t>
      void get_ds(
	lgr::opts_init_t<real_t> *arg
      )
      {
	throw std::runtime_error("dry_sizes does not feature a getter yet - TODO");
      }

      template <typename real_t>
      void get_dd(
	lgr::opts_init_t<real_t> *arg
      )
      {
	throw std::runtime_error("dry_distros does not feature a getter yet - TODO");
      }

      template <typename real_t>
      void get_sdd(
	lgr::opts_init_t<real_t> *arg
      )
      {
	throw std::runtime_error("source_dry_distros does not feature a getter yet - TODO");
      }

      template <typename real_t>
      void set_kp(
	lgr::opts_init_t<real_t> *arg,
	const bp::numeric::array &vec
      )
      {
        sanity_checks(vec);
	for (int i = 0; i < len(vec); ++i)
	  arg->kernel_parameters.push_back(bp::extract<real_t>(vec[i]));
      }

      template <typename real_t>
      bp::numeric::array get_kp(
	lgr::opts_init_t<real_t> *arg
      )
      {
	throw std::runtime_error("kernel_paramteres does not feature a getter yet - TODO");
      }
    };
  };
};
