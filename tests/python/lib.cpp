// written following:
// - http://www.boost.org/doc/libs/1_55_0/libs/python/doc/tutorial/doc/html/python/exposing.html
// - http://isolation-nation.blogspot.com/2008/09/packages-in-python-extension-modules.html

#include <boost/python.hpp> // TODO: less general to reduce compilation time
#include <boost/python/stl_iterator.hpp>

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_columnwise.hpp>

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

#include <libcloudph++/lgrngn/options.hpp>

using real_t = double;
namespace bp = boost::python;

struct iterable 
{
  bp::object &o;
  iterable(bp::object &o) : o(o) {}
  bp::stl_input_iterator<real_t> begin() const { return bp::stl_input_iterator<real_t>(o); }
  bp::stl_input_iterator<real_t> end() const { return bp::stl_input_iterator<real_t>(); }
};

namespace blk_1m
{
  void adj_cellwise(
    const libcloudphxx::blk_1m::opts_t<real_t>& opts,
    bp::object &rhod,
    bp::object &th,
    bp::object &rv,
    bp::object &rc,
    bp::object &rr,
    const real_t &dt
  )
  {
    iterable 
      th_c(th),
      rv_c(rv),
      rc_c(rc),
      rr_c(rr);
    libcloudphxx::blk_1m::adj_cellwise<real_t, iterable>(
      opts, 
      iterable(rhod), // since it is const, it may be a temporary object
      th_c, 
      rv_c, 
      rc_c, 
      rr_c, 
      dt
    );
  }

  void rhs_cellwise(
    const libcloudphxx::blk_1m::opts_t<real_t> &opts,
    bp::object &dot_rc,
    bp::object &dot_rr,
    bp::object &rc,
    bp::object &rr
  ) 
  {
    iterable 
      dot_rc_c(dot_rc), 
      dot_rr_c(dot_rr);
    libcloudphxx::blk_1m::rhs_cellwise(
      opts,
      dot_rc_c,
      dot_rr_c,
      iterable(rc),
      iterable(rr)
    );
  } 

  void rhs_columnwise(
    const libcloudphxx::blk_1m::opts_t<real_t> &opts,
    bp::object &dot_rr,
    bp::object &rhod,
    bp::object &rr,
    const real_t &dz
  ) {
    iterable 
      dot_rr_c(dot_rr);
// TODO: fails to compile (iterators)
//    libcloudphxx::blk_1m::rhs_columnwise<real_t, iterable>(
//      opts,
//      dot_rr_c,
//      iterable(rhod),
//      iterable(rr),
//      dz
//    );
  } 
};

namespace blk_2m
{
  void rhs_cellwise(
    const libcloudphxx::blk_2m::opts_t<real_t> &opts,
    bp::object &dot_th,
    bp::object &dot_rv,
    bp::object &dot_rc,
    bp::object &dot_nc,
    bp::object &dot_rr,
    bp::object &dot_nr,
    bp::object &rhod,
    bp::object &th,
    bp::object &rv,
    bp::object &rc,
    bp::object &nc,
    bp::object &rr,
    bp::object &nr,
    const real_t dt
  ) 
  {
    iterable 
      dot_th_c(dot_rc), 
      dot_rv_c(dot_rc), 
      dot_rc_c(dot_rc), 
      dot_nc_c(dot_rc), 
      dot_rr_c(dot_rc), 
      dot_nr_c(dot_rr);
// TODO: fails to compile (min() not declared for iterable)
/*
    libcloudphxx::blk_2m::rhs_cellwise(
      opts,
      dot_th_c,
      dot_rv_c,
      dot_rc_c,
      dot_nc_c,
      dot_rr_c,
      dot_nr_c,
      iterable(rhod),
      iterable(th),
      iterable(rv),
      iterable(rc),
      iterable(nc),
      iterable(rr),
      iterable(nr),
      dt
    );
*/
  } 

  void rhs_columnwise(
    const libcloudphxx::blk_2m::opts_t<real_t> &opts,
    bp::object &dot_rr,
    bp::object &dot_nr,
    bp::object &rhod,
    bp::object &rr,
    bp::object &nr,
    const real_t &dt,
    const real_t &dz
  ) {
    iterable 
      dot_rr_c(dot_rr), 
      dot_nr_c(dot_nr);
// TODO: fails to compile (iterators)
//    libcloudphxx::blk_2m::rhs_columnwise<real_t, iterable>(
//      opts,
//      dot_rr_c,
//      dot_nr_c,
//      iterable(rhod),
//      iterable(rr),
//      iterable(nr),
//      dt,
//      dz
//    );
  } 
};

namespace lgrngn
{
};

BOOST_PYTHON_MODULE(libcloudphxx)
{
  //bp::numeric::array::set_module_and_type("numpy", "ndarray");

  // specify that this module is actually a package
  bp::object package = bp::scope();
  package.attr("__path__") = "libcloudphxx";

  // blk_1m stuff
  {
    std::string nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".blk_1m");
    bp::object nested_module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
    bp::scope().attr("blk_1m") = nested_module;
    bp::scope parent = nested_module;
    bp::class_<libcloudphxx::blk_1m::opts_t<real_t>>("opts_t");
    bp::def("adj_cellwise", blk_1m::adj_cellwise);
    bp::def("rhs_cellwise", blk_1m::rhs_cellwise); 
    bp::def("rhs_columnwise", blk_1m::rhs_columnwise); 
  }

  // blk_2m stuff
  {
    std::string nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".blk_2m");
    bp::object nested_module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
    bp::scope().attr("blk_2m") = nested_module;
    bp::scope parent = nested_module;
    bp::class_<libcloudphxx::blk_2m::opts_t<real_t>>("opts_t");
    bp::def("rhs_cellwise", blk_2m::rhs_cellwise);
    bp::def("rhs_columnwise", blk_2m::rhs_columnwise);
  } 

  // lgrngn stuff
  {
    std::string nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".lgrngn");
    bp::object nested_module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
    bp::scope().attr("lgrngn") = nested_module;
    bp::scope parent = nested_module;
    bp::class_<libcloudphxx::lgrngn::opts_t<real_t>>("opts_t");
  }
}
