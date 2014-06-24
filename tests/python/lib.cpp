// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

#include <boost/python.hpp> 
// written following:
// - http://www.boost.org/doc/libs/1_55_0/libs/python/doc/tutorial/doc/html/python/exposing.html
// - http://isolation-nation.blogspot.com/2008/09/packages-in-python-extension-modules.html
#include <boost/assign/ptr_map_inserter.hpp>  // for 'ptr_map_insert()'

#if defined(BZ_THREADSAFE)
#  error please unset BZ_THREADSAFE
#endif
#include <blitz/array.h>

// all libcloudph++'s includes which potentially include <cassert> or <assert.h>
// needed here as assert.h redefines assert() every time it is included
#include <libcloudph++/blk_1m/extincl.hpp>
#include <libcloudph++/blk_2m/extincl.hpp>
#include <libcloudph++/lgrngn/extincl.hpp> // TODO: still missing some includes as asserts give aborts :(

// turning asserts into exceptions
#undef assert
#define assert_str(s) #s
// TODO: for some reason it's not working everywhere 
//#if defined(__GNUC_PREREQ) && __GNUC_PREREQ(2,6)
//#  define assert_fun __PRETTY_FUNCTION__
//#else
#  define assert_fun __func__
//#endif
#define assert(cond) { \
  if (!(cond)) throw std::runtime_error( \
    __FILE__ ":" + std::to_string(__LINE__) + ": " + assert_fun + ": " \
    "Assertion `" assert_str(cond) "' failed." \
  ); \
}

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_columnwise.hpp>

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

#include <libcloudph++/lgrngn/factory.hpp>

namespace bp = boost::python;
namespace b1m = libcloudphxx::blk_1m;
namespace b2m = libcloudphxx::blk_2m;
namespace lgr = libcloudphxx::lgrngn;
namespace cmn = libcloudphxx::common;

using real_t = double;
using arr_t = blitz::Array<real_t, 1>;
using py_ptr_t = long; // TODO: acquire it using some decltype()

// TODO: it does not when swithing from "numpy.array" to "numpy.ndarray"!
inline void sanity_checks(const bp::numeric::array &arg)
{
  // assuring double precision
  if (std::string(bp::extract<std::string>(arg.attr("dtype").attr("name"))) != "float64")
    throw std::runtime_error("dtype=float64 required for all passed arrays");

  // assuring contiguous layout
  if (!bp::extract<bool>(arg.attr("flags").attr("c_contiguous")))
    throw std::runtime_error("contiguous memory layout required");
}

inline arr_t np2bz(const bp::numeric::array &arg)
{
  sanity_checks(arg);

  // wrapping the data into a Blitz++ array to get STL-container-like functionality
  return arr_t(
    // pointer to the data
    reinterpret_cast<real_t*>(
      (py_ptr_t)bp::extract<py_ptr_t>(arg.attr("ctypes").attr("data")) 
    ), 
    // length of the array (regardless of the original dimensionality, we do 1D)
    blitz::shape(bp::extract<long>(arg.attr("size"))), 
    // ensure Blitz++ does not try to free the memory when done
    blitz::neverDeleteData
  );
}

const ptrdiff_t one = 1;

inline lgr::arrinfo_t<real_t> np2ai(const bp::numeric::array &arg)
{
  sanity_checks(arg);
  
  return lgr::arrinfo_t<real_t>(
    reinterpret_cast<real_t*>(
      (py_ptr_t)bp::extract<py_ptr_t>(arg.attr("ctypes").attr("data"))
    ),
    &one // TODO: parcel assumption hardcoded
  );
}

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
    return bp::object(bp::handle<>(PyBuffer_FromMemory(arg->outbuf(), len * sizeof(real_t))));
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
    const bp::numeric::array &rv
  )
  {
    lgr::arrinfo_t<real_t>
      np2ai_th(np2ai(th)),
      np2ai_rv(np2ai(rv));
    arg->step_sync(
      opts, 
      np2ai_th,
      np2ai_rv
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

namespace common
{
  template <typename real_t>
  real_t th_dry2std(
    const real_t &th_dry,
    const real_t &r
  )
  {
    return cmn::theta_dry::dry2std(th_dry * si::kelvins, r * si::dimensionless()) / si::kelvins;
  }

  template <typename real_t>
  real_t th_std2dry(
    const real_t &th_std,
    const real_t &r
  )
  {
    return cmn::theta_dry::std2dry(th_std * si::kelvins, r * si::dimensionless()) / si::kelvins;
  }

  template <typename real_t>
  real_t p_vs(const real_t &T)
  {
    return cmn::const_cp::p_vs(T * si::kelvins) / si::pascals;
  }

  template <typename real_t>
  real_t T(const real_t &th, const real_t &rhod)
  {
    return cmn::theta_dry::T(th * si::kelvins, rhod  * si::kilograms / si::cubic_metres) / si::kelvins;
  }

  template <typename real_t>
  real_t p(const real_t &rhod, const real_t &r, const real_t &T)
  {
    return cmn::theta_dry::p(rhod  * si::kilograms / si::cubic_metres, r * si::kilograms / si::kilograms, T * si::kelvins) / si::pascals;
  }
};


BOOST_PYTHON_MODULE(libcloudphxx)
{
  bp::numeric::array::set_module_and_type("numpy", "ndarray");

  // specify that this module is actually a package
  bp::object package = bp::scope();
  package.attr("__path__") = "libcloudphxx";

  // common stuff
  {
    std::string nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".common");
    bp::object nested_module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
    bp::scope().attr("common") = nested_module;
    bp::scope parent = nested_module;

    bp::scope().attr("R_d") = (real_t) (cmn::moist_air::R_d<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("c_pd") = (real_t) (cmn::moist_air::c_pd<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("g") = (real_t) (cmn::earth::g<real_t>() / si::metres * si::seconds * si::seconds);
    bp::scope().attr("p_1000") = (real_t) (cmn::theta_std::p_1000<real_t>() / si::pascals);
    bp::scope().attr("eps") = (real_t) (cmn::moist_air::eps<real_t>());
    // TODO: how to make the above constant?

    bp::def("th_dry2std", &common::th_dry2std<real_t>);
    bp::def("th_std2dry", &common::th_std2dry<real_t>);
    bp::def("p_vs", &common::p_vs<real_t>);
    bp::def("T", &common::T<real_t>);
    bp::def("p", &common::p<real_t>);
  }

  // blk_1m stuff
  {
    std::string nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".blk_1m");
    bp::object nested_module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
    bp::scope().attr("blk_1m") = nested_module;
    bp::scope parent = nested_module;
    bp::class_<b1m::opts_t<real_t>>("opts_t")
      .def_readwrite("cond", &b1m::opts_t<real_t>::cond)
      .def_readwrite("cevp", &b1m::opts_t<real_t>::cevp)
      .def_readwrite("revp", &b1m::opts_t<real_t>::revp)
      .def_readwrite("conv", &b1m::opts_t<real_t>::conv)
      .def_readwrite("accr", &b1m::opts_t<real_t>::accr)
      .def_readwrite("sedi", &b1m::opts_t<real_t>::sedi)
      .def_readwrite("r_c0", &b1m::opts_t<real_t>::r_c0)
      .def_readwrite("r_eps", &b1m::opts_t<real_t>::r_eps)
      ;
    bp::def("adj_cellwise", blk_1m::adj_cellwise);
    bp::def("rhs_cellwise", blk_1m::rhs_cellwise); 
    bp::def("rhs_columnwise", blk_1m::rhs_columnwise); // TODO: handle the returned flux
  }

  // blk_2m stuff
  {
    std::string nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".blk_2m");
    bp::object nested_module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
    bp::scope().attr("blk_2m") = nested_module;
    bp::scope parent = nested_module;
    bp::class_<b2m::opts_t<real_t>>("opts_t")
      .def_readwrite("acti", &b2m::opts_t<real_t>::acti)
      .def_readwrite("cond", &b2m::opts_t<real_t>::cond)
      .def_readwrite("acnv", &b2m::opts_t<real_t>::acnv)
      .def_readwrite("accr", &b2m::opts_t<real_t>::accr)
      .def_readwrite("sedi", &b2m::opts_t<real_t>::sedi)
      .def_readwrite("RH_max", &b2m::opts_t<real_t>::RH_max)
      .add_property("dry_distros", &blk_2m::get_dd, &blk_2m::set_dd)
    ;
    bp::def("rhs_cellwise", blk_2m::rhs_cellwise);
    bp::def("rhs_columnwise", blk_2m::rhs_columnwise); // TODO: handle the returned flux
  } 

  // lgrngn stuff
  {
    std::string nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".lgrngn");
    bp::object nested_module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
    bp::scope().attr("lgrngn") = nested_module;
    bp::scope parent = nested_module;
    // enums
    bp::enum_<lgr::backend_t>("backend_t")
      .value("serial", lgr::serial)
      .value("OpenMP", lgr::OpenMP)
      .value("CUDA",   lgr::CUDA);
    bp::enum_<lgr::kernel_t>("kernel_t") 
      .value("geometric", lgr::geometric);

    bp::enum_<lgr::chem_gas>("chem_gas")
      .value("gSO2",  lgr::gSO2)
      .value("gO3",   lgr::gO3)
      .value("gH2O2", lgr::gH2O2);

    bp::enum_<lgr::chem_aq>("chem_aq")
      .value("H",    lgr::H)
      .value("OH",   lgr::OH)
      .value("SO2",  lgr::SO2)
      .value("O3",   lgr::O3)
      .value("H2O2", lgr::H2O2)
      .value("HSO3", lgr::HSO3)
      .value("SO3",  lgr::SO3)
      .value("S_VI", lgr::S_VI)
      .value("HSO4", lgr::HSO4)
      .value("SO4",  lgr::SO4);
    // classes
    bp::class_<lgr::opts_t<real_t>>("opts_t")
      .def_readwrite("adve", &lgr::opts_t<real_t>::adve)
      .def_readwrite("sedi", &lgr::opts_t<real_t>::sedi)
      .def_readwrite("cond", &lgr::opts_t<real_t>::cond)
      .def_readwrite("coal", &lgr::opts_t<real_t>::coal)
      .def_readwrite("chem", &lgr::opts_t<real_t>::chem)
      .def_readwrite("RH_max", &lgr::opts_t<real_t>::RH_max)
      .def_readwrite("sstp_cond", &lgr::opts_t<real_t>::sstp_cond)
      .def_readwrite("sstp_coal", &lgr::opts_t<real_t>::sstp_coal)
      .def_readwrite("sstp_chem", &lgr::opts_t<real_t>::sstp_chem)
      .add_property("chem_gas", &lgrngn::get_cg, &lgrngn::set_cg)
    ;
    bp::class_<lgr::opts_init_t<real_t>>("opts_init_t")
      .add_property("dry_distros", &lgrngn::get_dd, &lgrngn::set_dd)
      .def_readwrite("nx", &lgr::opts_init_t<real_t>::nx)
      .def_readwrite("ny", &lgr::opts_init_t<real_t>::ny)
      .def_readwrite("nz", &lgr::opts_init_t<real_t>::nz)
      .def_readwrite("dx", &lgr::opts_init_t<real_t>::dx)
      .def_readwrite("dy", &lgr::opts_init_t<real_t>::dy)
      .def_readwrite("dz", &lgr::opts_init_t<real_t>::dz)
      .def_readwrite("x0", &lgr::opts_init_t<real_t>::x0)
      .def_readwrite("x1", &lgr::opts_init_t<real_t>::x1)
      .def_readwrite("y0", &lgr::opts_init_t<real_t>::y0)
      .def_readwrite("y1", &lgr::opts_init_t<real_t>::y1)
      .def_readwrite("z0", &lgr::opts_init_t<real_t>::z0)
      .def_readwrite("z1", &lgr::opts_init_t<real_t>::z1)
      .def_readwrite("dt", &lgr::opts_init_t<real_t>::dt)
      .def_readwrite("kernel", &lgr::opts_init_t<real_t>::kernel)
      .def_readwrite("sd_conc_mean", &lgr::opts_init_t<real_t>::sd_conc_mean)
      .def_readwrite("chem_rho", &lgr::opts_init_t<real_t>::chem_rho)
    ;
    bp::class_<lgr::particles_proto_t<real_t>/*, boost::noncopyable*/>("particles_proto_t")
      .def("init",         &lgrngn::init)
      .def("step_sync",    &lgrngn::step_sync)
      .def("step_async",   &lgr::particles_proto_t<real_t>::step_async)
      .def("diag_sd_conc", &lgr::particles_proto_t<real_t>::diag_sd_conc)
      .def("diag_dry_rng", &lgr::particles_proto_t<real_t>::diag_dry_rng)
      .def("diag_wet_rng", &lgr::particles_proto_t<real_t>::diag_wet_rng)
      .def("diag_dry_mom", &lgr::particles_proto_t<real_t>::diag_dry_mom)
      .def("diag_wet_mom", &lgr::particles_proto_t<real_t>::diag_wet_mom)
      .def("diag_chem",    &lgr::particles_proto_t<real_t>::diag_chem)
      .def("outbuf",       &lgrngn::outbuf)
    ;
    // functions
    bp::def("factory", lgrngn::factory, bp::return_value_policy<bp::manage_new_object>());
  }
}
