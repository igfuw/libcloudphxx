// written following:
// - http://www.boost.org/doc/libs/1_55_0/libs/python/doc/tutorial/doc/html/python/exposing.html
// - http://isolation-nation.blogspot.com/2008/09/packages-in-python-extension-modules.html

#include <boost/python.hpp> // TODO: less general to reduce compilation time
#include <boost/python/stl_iterator.hpp>

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>

using real_t = double;
namespace bp = boost::python;

struct iterable 
{
  bp::object &o;
  iterable(bp::object &o) : o(o) {}
  bp::stl_input_iterator<real_t> begin() const { return bp::stl_input_iterator<real_t>(o); }
  bp::stl_input_iterator<real_t> end() const { return bp::stl_input_iterator<real_t>(); }
};

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
  iterable th_c(th), rv_c(rv), rc_c(rc), rr_c(rr);
  libcloudphxx::blk_1m::adj_cellwise<real_t, iterable>(
    opts, 
    iterable(rhod), 
    th_c, 
    rv_c, 
    rc_c, 
    rr_c, 
    dt
  );
}

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
    bp::def("adj_cellwise", adj_cellwise);
  }

  // blk_2m stuff
  {} // TODO

  // lgrngn stuff
  {} // TODO
}
