// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw
//
// written with a lot of help from:
// - http://www.boost.org/doc/libs/1_55_0/libs/python/doc/tutorial/doc/html/python/exposing.html
// - http://isolation-nation.blogspot.com/2008/09/packages-in-python-extension-modules.html

#include <libcloudph++/git_revision.h>

#include "util.hpp"
#include "blk_1m.hpp"
#include "blk_2m.hpp"
#include "lgrngn.hpp"
#include "common.hpp"

BOOST_PYTHON_MODULE(libcloudphxx)
{
  namespace bp = boost::python;
  using namespace libcloudphxx::python;

  using real_t = double;
  using arr_t = blitz::Array<real_t, 1>;

  bp::numeric::array::set_module_and_type("numpy", "ndarray");

  // specify that this module is actually a package
  bp::object package = bp::scope();
  package.attr("__path__") = "libcloudphxx";

  // exposing git revision id
  package.attr("git_revision") = GIT_REVISION;

  // common stuff
  {
    std::string nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".common");
    bp::object nested_module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
    bp::scope().attr("common") = nested_module;
    bp::scope parent = nested_module;

    bp::scope().attr("R_d") = (real_t) (cmn::moist_air::R_d<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("R_v") = (real_t) (cmn::moist_air::R_v<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("c_pd") = (real_t) (cmn::moist_air::c_pd<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("c_pv") = (real_t) (cmn::moist_air::c_pv<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("g") = (real_t) (cmn::earth::g<real_t>() / si::metres * si::seconds * si::seconds);
    bp::scope().attr("p_1000") = (real_t) (cmn::theta_std::p_1000<real_t>() / si::pascals);
    bp::scope().attr("eps") = (real_t) (cmn::moist_air::eps<real_t>());
    bp::scope().attr("rho_w") = (real_t) (cmn::moist_air::rho_w<real_t>() * si::cubic_metres / si::kilograms);
    // TODO: how to make the above constant?

    bp::def("th_dry2std", &common::th_dry2std<real_t>);
    bp::def("th_std2dry", &common::th_std2dry<real_t>);
    bp::def("p_vs", &common::p_vs<real_t>);
    bp::def("T", &common::T<real_t>);
    bp::def("p", &common::p<real_t>);
    bp::def("rw3_cr", &common::rw3_cr<real_t>);
    bp::def("S_cr", &common::S_cr<real_t>);
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
    bp::def("adj_cellwise", blk_1m::adj_cellwise<arr_t>);
    bp::def("rhs_cellwise", blk_1m::rhs_cellwise<arr_t>); 
    bp::def("rhs_columnwise", blk_1m::rhs_columnwise<arr_t>); // TODO: handle the returned flux
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
      .add_property("dry_distros", &blk_2m::get_dd<real_t>, &blk_2m::set_dd<real_t>)
    ;
    bp::def("rhs_cellwise", blk_2m::rhs_cellwise<arr_t>);
    bp::def("rhs_columnwise", blk_2m::rhs_columnwise<arr_t>); // TODO: handle the returned flux
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
      .value("geometric", lgr::geometric)
      .value("golovin", lgr::golovin);

    bp::enum_<lgr::chem_species_t>("chem_species_t")
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
      .add_property("chem_gas", &lgrngn::get_cg<real_t>, &lgrngn::set_cg<real_t>)
    ;
    bp::class_<lgr::opts_init_t<real_t>>("opts_init_t")
      .add_property("dry_distros", &lgrngn::get_dd<real_t>, &lgrngn::set_dd<real_t>)
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
      .def_readwrite("chem_switch", &lgr::opts_init_t<real_t>::chem_switch)
      .def_readwrite("sstp_cond", &lgr::opts_init_t<real_t>::sstp_cond)
      .def_readwrite("sstp_coal", &lgr::opts_init_t<real_t>::sstp_coal)
      .def_readwrite("sstp_chem", &lgr::opts_init_t<real_t>::sstp_chem)
      .def_readwrite("kernel", &lgr::opts_init_t<real_t>::kernel)
      .def_readwrite("sd_conc_mean", &lgr::opts_init_t<real_t>::sd_conc_mean)
      .def_readwrite("chem_rho", &lgr::opts_init_t<real_t>::chem_rho)
      .def_readwrite("RH_max", &lgr::opts_init_t<real_t>::RH_max)
      .add_property("kernel_parameters", &lgrngn::get_kp<real_t>, &lgrngn::set_kp<real_t>)
    ;
    bp::class_<lgr::particles_proto_t<real_t>/*, boost::noncopyable*/>("particles_proto_t")
      .add_property("opts_init", &lgrngn::get_oi<real_t>)
      .def("init",         &lgrngn::init_3arg<real_t>)
      .def("init",         &lgrngn::init_5arg<real_t>)
      .def("step_sync",    &lgrngn::step_sync<real_t>)
      .def("step_async",   &lgr::particles_proto_t<real_t>::step_async)
      .def("diag_sd_conc", &lgr::particles_proto_t<real_t>::diag_sd_conc)
      .def("diag_all",     &lgr::particles_proto_t<real_t>::diag_all)
      .def("diag_rw_ge_rc",&lgr::particles_proto_t<real_t>::diag_rw_ge_rc)
      .def("diag_dry_rng", &lgr::particles_proto_t<real_t>::diag_dry_rng)
      .def("diag_wet_rng", &lgr::particles_proto_t<real_t>::diag_wet_rng)
      .def("diag_dry_mom", &lgr::particles_proto_t<real_t>::diag_dry_mom)
      .def("diag_wet_mom", &lgr::particles_proto_t<real_t>::diag_wet_mom)
      .def("diag_wet_mass_dens", &lgr::particles_proto_t<real_t>::diag_wet_mass_dens)
      .def("diag_chem",    &lgr::particles_proto_t<real_t>::diag_chem)
      .def("outbuf",       &lgrngn::outbuf<real_t>)
    ;
    // functions
    bp::def("factory", lgrngn::factory<real_t>, bp::return_value_policy<bp::manage_new_object>());
  }
}
