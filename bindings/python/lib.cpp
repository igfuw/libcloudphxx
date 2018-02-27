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


#ifdef BPNUMERIC
  #define BP_ARR_FROM_BP_OBJ  bp_array(bp::object())
#elif defined BPNUMPY
  #define BP_ARR_FROM_BP_OBJ  bp_array(bp::numpy::array(bp::object()))
#endif

BOOST_PYTHON_MODULE(libcloudphxx)
{
  namespace bp = boost::python;
  using namespace libcloudphxx::python;

  using real_t = double;
  using arr_t = blitz::Array<real_t, 1>;

#ifdef BPNUMERIC
  bp_array::set_module_and_type("numpy", "ndarray");
#elif defined BPNUMPY
  Py_Initialize();
  bp::numpy::initialize();
#endif

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

    bp::scope().attr("R") = (real_t) (cmn::moist_air::kaBoNA<real_t>() / si::joules * si::kelvins * si::moles);
    bp::scope().attr("R_d") = (real_t) (cmn::moist_air::R_d<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("R_v") = (real_t) (cmn::moist_air::R_v<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("c_pd") = (real_t) (cmn::moist_air::c_pd<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("c_pv") = (real_t) (cmn::moist_air::c_pv<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("c_pw") = (real_t) (cmn::moist_air::c_pw<real_t>() / si::joules * si::kilograms * si::kelvins);
    bp::scope().attr("g") = (real_t) (cmn::earth::g<real_t>() / si::metres * si::seconds * si::seconds);
    bp::scope().attr("p_1000") = (real_t) (cmn::theta_std::p_1000<real_t>() / si::pascals);
    bp::scope().attr("eps") = (real_t) (cmn::moist_air::eps<real_t>());
    bp::scope().attr("rho_stp") = (real_t) (cmn::earth::rho_stp<real_t>() * si::cubic_metres / si::kilograms);
    bp::scope().attr("rho_w") = (real_t) (cmn::moist_air::rho_w<real_t>() * si::cubic_metres / si::kilograms);
    //molar mass of trace gases
    bp::scope().attr("M_SO2")   = (real_t) (cmn::molar_mass::M_SO2<real_t>()   * si::moles / si::kilograms);
    bp::scope().attr("M_H2O2")  = (real_t) (cmn::molar_mass::M_H2O2<real_t>()  * si::moles / si::kilograms);
    bp::scope().attr("M_O3")    = (real_t) (cmn::molar_mass::M_O3<real_t>()    * si::moles / si::kilograms);
    bp::scope().attr("M_NH3")   = (real_t) (cmn::molar_mass::M_NH3<real_t>()   * si::moles / si::kilograms);
    bp::scope().attr("M_HNO3")  = (real_t) (cmn::molar_mass::M_HNO3<real_t>()  * si::moles / si::kilograms);
    bp::scope().attr("M_CO2")   = (real_t) (cmn::molar_mass::M_CO2<real_t>()    * si::moles / si::kilograms);
    //H+ and OH-
    bp::scope().attr("M_H")     = (real_t) (cmn::molar_mass::M_H<real_t>()     * si::moles / si::kilograms);
    bp::scope().attr("M_OH")    = (real_t) (cmn::molar_mass::M_OH<real_t>()    * si::moles / si::kilograms);
    //SO2 * H20
    bp::scope().attr("M_SO2_H2O")   = (real_t) (cmn::molar_mass::M_SO2_H2O<real_t>()   * si::moles / si::kilograms);
    bp::scope().attr("M_HSO3")  = (real_t) (cmn::molar_mass::M_HSO3<real_t>()  * si::moles / si::kilograms);
    bp::scope().attr("M_SO3")   = (real_t) (cmn::molar_mass::M_SO3<real_t>()   * si::moles / si::kilograms);
    //HSO4
    bp::scope().attr("M_H2SO4") = (real_t) (cmn::molar_mass::M_H2SO4<real_t>() * si::moles / si::kilograms);
    bp::scope().attr("M_HSO4")  = (real_t) (cmn::molar_mass::M_HSO4<real_t>()  * si::moles / si::kilograms);
    bp::scope().attr("M_SO4")   = (real_t) (cmn::molar_mass::M_SO4<real_t>()   * si::moles / si::kilograms);
    //CO2 * H20
    bp::scope().attr("M_CO2_H2O") = (real_t) (cmn::molar_mass::M_CO2_H2O<real_t>()   * si::moles / si::kilograms);
    bp::scope().attr("M_HCO3")    = (real_t) (cmn::molar_mass::M_HCO3<real_t>()  * si::moles / si::kilograms);
    bp::scope().attr("M_CO3")     = (real_t) (cmn::molar_mass::M_CO3<real_t>()   * si::moles / si::kilograms);
    //NH3 * H20
    bp::scope().attr("M_NH3_H2O") = (real_t) (cmn::molar_mass::M_NH3_H2O<real_t>()   * si::moles / si::kilograms);
    bp::scope().attr("M_NH4")     = (real_t) (cmn::molar_mass::M_NH4<real_t>()  * si::moles / si::kilograms);
    //HNO3(aq)
    bp::scope().attr("M_NO3")     = (real_t) (cmn::molar_mass::M_NO3<real_t>()   * si::moles / si::kilograms);
 
    //Henry const
    bp::scope().attr("H_SO2")  = (real_t) (cmn::henry::H_SO2<real_t>()  * si::cubic_metres * si::pascals / si::moles);
    bp::scope().attr("H_H2O2") = (real_t) (cmn::henry::H_H2O2<real_t>() * si::cubic_metres * si::pascals / si::moles);
    bp::scope().attr("H_O3")   = (real_t) (cmn::henry::H_O3<real_t>()   * si::cubic_metres * si::pascals / si::moles);
    bp::scope().attr("H_HNO3") = (real_t) (cmn::henry::H_HNO3<real_t>() * si::cubic_metres * si::pascals / si::moles);
    bp::scope().attr("H_NH3")  = (real_t) (cmn::henry::H_NH3<real_t>()  * si::cubic_metres * si::pascals / si::moles);
    bp::scope().attr("H_CO2")  = (real_t) (cmn::henry::H_CO2<real_t>()  * si::cubic_metres * si::pascals / si::moles);
    //Henry dep on temperature
    bp::scope().attr("dHR_SO2")  = (real_t) (cmn::henry::dHR_SO2<real_t>()  / si::kelvins);
    bp::scope().attr("dHR_H2O2") = (real_t) (cmn::henry::dHR_H2O2<real_t>() / si::kelvins);
    bp::scope().attr("dHR_O3")   = (real_t) (cmn::henry::dHR_O3<real_t>()   / si::kelvins);
    bp::scope().attr("dHR_HNO3") = (real_t) (cmn::henry::dHR_HNO3<real_t>() / si::kelvins);
    bp::scope().attr("dHR_NH3")  = (real_t) (cmn::henry::dHR_NH3<real_t>()  / si::kelvins);
    bp::scope().attr("dHR_CO2")  = (real_t) (cmn::henry::dHR_CO2<real_t>()  / si::kelvins);
 
    //dissoc const
    bp::scope().attr("K_H2O")  = (real_t) (cmn::dissoc::K_H2O<real_t>()  / si::moles / si::moles * si::cubic_metres * si::cubic_metres);
    bp::scope().attr("K_SO2")  = (real_t) (cmn::dissoc::K_SO2<real_t>()  / si::moles * si::cubic_metres);
    bp::scope().attr("K_HSO3") = (real_t) (cmn::dissoc::K_HSO3<real_t>() / si::moles * si::cubic_metres);
    bp::scope().attr("K_HSO4") = (real_t) (cmn::dissoc::K_HSO4<real_t>() / si::moles * si::cubic_metres);
    bp::scope().attr("K_CO2")  = (real_t) (cmn::dissoc::K_CO2<real_t>()  / si::moles * si::cubic_metres);
    bp::scope().attr("K_HCO3") = (real_t) (cmn::dissoc::K_HCO3<real_t>() / si::moles * si::cubic_metres);
    bp::scope().attr("K_NH3")  = (real_t) (cmn::dissoc::K_NH3<real_t>()  / si::moles * si::cubic_metres);
    bp::scope().attr("K_HNO3") = (real_t) (cmn::dissoc::K_HNO3<real_t>() / si::moles * si::cubic_metres);
    //dissoc const dep on temperature
    bp::scope().attr("dKR_SO2")  = (real_t) (cmn::dissoc::dKR_SO2<real_t>()  / si::kelvins);
    bp::scope().attr("dKR_HSO3") = (real_t) (cmn::dissoc::dKR_HSO3<real_t>() / si::kelvins);
    bp::scope().attr("dKR_HSO4") = (real_t) (cmn::dissoc::dKR_HSO4<real_t>() / si::kelvins);
    bp::scope().attr("dKR_HNO3") = (real_t) (cmn::dissoc::dKR_HNO3<real_t>() / si::kelvins);
    bp::scope().attr("dKR_NH3")  = (real_t) (cmn::dissoc::dKR_NH3<real_t>()  / si::kelvins);
    bp::scope().attr("dKR_CO2")  = (real_t) (cmn::dissoc::dKR_CO2<real_t>()  / si::kelvins);
    bp::scope().attr("dKR_HCO3") = (real_t) (cmn::dissoc::dKR_HCO3<real_t>() / si::kelvins);
    // TODO: how to make the above constant?

    bp::def("th_dry2std", &common::th_dry2std<real_t>);
    bp::def("th_std2dry", &common::th_std2dry<real_t>);
    bp::def("p_vs", &common::p_vs<real_t>);
    bp::def("l_v", &common::l_v<real_t>);
    bp::def("T", &common::T<real_t>);
    bp::def("p", &common::p<real_t>);
    bp::def("rw3_cr", &common::rw3_cr<real_t>);
    bp::def("S_cr", &common::S_cr<real_t>);
    bp::def("p_hydro", &common::p_hydro<real_t>);
    bp::def("rhod", &common::rhod<real_t>);
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
    bp::enum_<lgr::kernel_t::kernel_t>("kernel_t") 
      .value("geometric", lgr::kernel_t::geometric)
      .value("golovin", lgr::kernel_t::golovin)
      .value("hall", lgr::kernel_t::hall)
      .value("hall_davis_no_waals", lgr::kernel_t::hall_davis_no_waals)
      .value("long", lgr::kernel_t::Long)
      .value("hall_pinsky_1000mb_grav", lgr::kernel_t::hall_pinsky_1000mb_grav)
      .value("onishi_hall_davis_no_waals", lgr::kernel_t::onishi_hall_davis_no_waals)
      .value("onishi_hall", lgr::kernel_t::onishi_hall)
      .value("vohl_davis_no_waals", lgr::kernel_t::vohl_davis_no_waals)
      .value("hall_pinsky_stratocumulus", lgr::kernel_t::hall_pinsky_stratocumulus)
      .value("hall_pinsky_cumulonimbus", lgr::kernel_t::hall_pinsky_cumulonimbus);
    bp::enum_<lgr::vt_t::vt_t>("vt_t") 
      .value("beard76", lgr::vt_t::beard76)
      .value("beard77", lgr::vt_t::beard77)
      .value("beard77fast", lgr::vt_t::beard77fast)
      .value("khvorostyanov_spherical", lgr::vt_t::khvorostyanov_spherical)
      .value("khvorostyanov_nonspherical", lgr::vt_t::khvorostyanov_nonspherical);
    bp::enum_<lgr::as_t::as_t>("as_t") 
      .value("implicit", lgr::as_t::implicit)
      .value("euler", lgr::as_t::euler)
      .value("pred_corr", lgr::as_t::pred_corr);

    bp::enum_<lgr::chem_species_t>("chem_species_t")
      .value("H",    lgr::H)
      .value("SO2",  lgr::SO2)
      .value("O3",   lgr::O3)
      .value("H2O2", lgr::H2O2)
      .value("CO2",  lgr::CO2)
      .value("NH3",  lgr::NH3)
      .value("HNO3", lgr::HNO3)
      .value("S_VI", lgr::S_VI);
    // classes
    bp::class_<lgr::opts_t<real_t>>("opts_t")
      .def_readwrite("adve", &lgr::opts_t<real_t>::adve)
      .def_readwrite("sedi", &lgr::opts_t<real_t>::sedi)
      .def_readwrite("cond", &lgr::opts_t<real_t>::cond)
      .def_readwrite("coal", &lgr::opts_t<real_t>::coal)
      .def_readwrite("src", &lgr::opts_t<real_t>::src)
      .def_readwrite("rcyc", &lgr::opts_t<real_t>::rcyc)
      .def_readwrite("chem_dsl", &lgr::opts_t<real_t>::chem_dsl)
      .def_readwrite("chem_dsc", &lgr::opts_t<real_t>::chem_dsc)
      .def_readwrite("chem_rct", &lgr::opts_t<real_t>::chem_rct)
      .def_readwrite("RH_max", &lgr::opts_t<real_t>::RH_max)
    ;
    bp::class_<lgr::opts_init_t<real_t>>("opts_init_t")
      .add_property("dry_distros", &lgrngn::get_dd<real_t>, &lgrngn::set_dd<real_t>)
      .add_property("src_dry_distros", &lgrngn::get_sdd<real_t>, &lgrngn::set_sdd<real_t>)
      .add_property("dry_sizes", &lgrngn::get_ds<real_t>, &lgrngn::set_ds<real_t>)
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
      .def_readwrite("dev_id", &lgr::opts_init_t<real_t>::dev_id)
      .def_readwrite("dev_count", &lgr::opts_init_t<real_t>::dev_count)
      .def_readwrite("src_z1", &lgr::opts_init_t<real_t>::src_z1)
      .def_readwrite("dt", &lgr::opts_init_t<real_t>::dt)
      .def_readwrite("chem_switch", &lgr::opts_init_t<real_t>::chem_switch)
      .def_readwrite("coal_switch", &lgr::opts_init_t<real_t>::coal_switch)
      .def_readwrite("sedi_switch", &lgr::opts_init_t<real_t>::sedi_switch)
      .def_readwrite("src_switch", &lgr::opts_init_t<real_t>::src_switch)
      .def_readwrite("exact_sstp_cond", &lgr::opts_init_t<real_t>::exact_sstp_cond)
      .def_readwrite("sstp_cond", &lgr::opts_init_t<real_t>::sstp_cond)
      .def_readwrite("sstp_coal", &lgr::opts_init_t<real_t>::sstp_coal)
      .def_readwrite("sstp_chem", &lgr::opts_init_t<real_t>::sstp_chem)
      .def_readwrite("supstp_src", &lgr::opts_init_t<real_t>::supstp_src)
      .def_readwrite("kernel", &lgr::opts_init_t<real_t>::kernel)
      .def_readwrite("adve_scheme", &lgr::opts_init_t<real_t>::adve_scheme)
      .def_readwrite("sd_conc", &lgr::opts_init_t<real_t>::sd_conc)
      .def_readwrite("sd_conc_large_tail", &lgr::opts_init_t<real_t>::sd_conc_large_tail)
      .def_readwrite("sd_const_multi", &lgr::opts_init_t<real_t>::sd_const_multi)
      .def_readwrite("sd_const_multi_dry_sizes", &lgr::opts_init_t<real_t>::sd_const_multi_dry_sizes)
      .def_readwrite("src_sd_conc", &lgr::opts_init_t<real_t>::src_sd_conc)
      .def_readwrite("n_sd_max", &lgr::opts_init_t<real_t>::n_sd_max)
      .def_readwrite("terminal_velocity", &lgr::opts_init_t<real_t>::terminal_velocity)
      .def_readwrite("chem_rho", &lgr::opts_init_t<real_t>::chem_rho)
      .def_readwrite("RH_max", &lgr::opts_init_t<real_t>::RH_max)
      .def_readwrite("rng_seed", &lgr::opts_init_t<real_t>::rng_seed)
      .add_property("kernel_parameters", &lgrngn::get_kp<real_t>, &lgrngn::set_kp<real_t>)
    ;
    bp::class_<lgr::particles_proto_t<real_t>/*, boost::noncopyable*/>("particles_proto_t")
      .add_property("opts_init", &lgrngn::get_oi<real_t>)
      .def("init",         &lgrngn::init<real_t>, (
        bp::arg("th")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("rv")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("rhod")= BP_ARR_FROM_BP_OBJ,
        bp::arg("Cx")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("Cy")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("Cz")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("ambient_chem") = bp::dict()
      ))
      .def("step_sync",    &lgrngn::step_sync<real_t>, (
        bp::arg("th")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("rv")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("rhod")= BP_ARR_FROM_BP_OBJ,
        bp::arg("Cx")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("Cy")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("Cz")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("RH")  = BP_ARR_FROM_BP_OBJ,
        bp::arg("ambient_chem") = bp::dict()
      ))
      .def("step_async",   &lgr::particles_proto_t<real_t>::step_async)
      .def("diag_sd_conc", &lgr::particles_proto_t<real_t>::diag_sd_conc)
      .def("diag_all",     &lgr::particles_proto_t<real_t>::diag_all)
      .def("diag_rw_ge_rc",&lgr::particles_proto_t<real_t>::diag_rw_ge_rc)
      .def("diag_RH_ge_Sc",&lgr::particles_proto_t<real_t>::diag_RH_ge_Sc)
      .def("diag_RH",&lgr::particles_proto_t<real_t>::diag_RH)
      .def("diag_vel_div",&lgr::particles_proto_t<real_t>::diag_vel_div)
      .def("diag_dry_rng", &lgr::particles_proto_t<real_t>::diag_dry_rng)
      .def("diag_wet_rng", &lgr::particles_proto_t<real_t>::diag_wet_rng)
      .def("diag_kappa_rng", &lgr::particles_proto_t<real_t>::diag_kappa_rng)
      .def("diag_dry_mom", &lgr::particles_proto_t<real_t>::diag_dry_mom)
      .def("diag_wet_mom", &lgr::particles_proto_t<real_t>::diag_wet_mom)
      .def("diag_kappa_mom",    &lgr::particles_proto_t<real_t>::diag_kappa_mom)
      .def("diag_wet_mass_dens", &lgr::particles_proto_t<real_t>::diag_wet_mass_dens)
      .def("diag_chem",    &lgr::particles_proto_t<real_t>::diag_chem)
      .def("diag_precip_rate",    &lgr::particles_proto_t<real_t>::diag_precip_rate)
      .def("diag_puddle",    &lgrngn::diag_puddle<real_t>)
      .def("outbuf",       &lgrngn::outbuf<real_t>)
    ;
    // functions
    bp::def("factory", lgrngn::factory<real_t>, bp::return_value_policy<bp::manage_new_object>());
  }
}
