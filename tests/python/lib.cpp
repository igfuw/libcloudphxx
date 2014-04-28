// written following:
// - http://www.boost.org/doc/libs/1_55_0/libs/python/doc/tutorial/doc/html/python/exposing.html
// - http://isolation-nation.blogspot.com/2008/09/packages-in-python-extension-modules.html

#include <boost/python.hpp> 

#if defined(BZ_THREADSAFE)
#  error please unset BZ_THREADSAFE
#endif
#include <blitz/array.h>

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_columnwise.hpp>

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

#include <libcloudph++/lgrngn/opts.hpp>

namespace bp = boost::python;

using real_t = double;
using arr_t = blitz::Array<real_t, 1>;

inline arr_t np2bz(bp::numeric::array &arg)
{
// TODO: asserts on shape?
  if (std::string(bp::extract<std::string>(arg.attr("dtype").attr("name"))) != "float64")
    throw std::runtime_error("dtype=float64 required for all passed arrays");

  return arr_t(
    reinterpret_cast<real_t*>(
      (long)bp::extract<long>(arg.attr("ctypes").attr("data"))
    ), 
    blitz::shape(bp::extract<long>(arg.attr("size"))), 
    blitz::neverDeleteData
  );
}

namespace blk_1m
{
  void adj_cellwise(
    const libcloudphxx::blk_1m::opts_t<real_t>& opts,
    bp::numeric::array &rhod,
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
    libcloudphxx::blk_1m::adj_cellwise(
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
    const libcloudphxx::blk_1m::opts_t<real_t> &opts,
    bp::numeric::array &dot_rc,
    bp::numeric::array &dot_rr,
    bp::numeric::array &rc,
    bp::numeric::array &rr
  ) 
  {
    arr_t
      np2bz_dot_rc(np2bz(dot_rc)), 
      np2bz_dot_rr(np2bz(dot_rr));
    libcloudphxx::blk_1m::rhs_cellwise(
      opts,
      np2bz_dot_rc,
      np2bz_dot_rr,
      np2bz(rc),
      np2bz(rr)
    );
  } 

/*
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
*/
};

namespace blk_2m
{
/*
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
*/
};

namespace lgrngn
{
};

BOOST_PYTHON_MODULE(libcloudphxx)
{
  bp::numeric::array::set_module_and_type("numpy", "ndarray");

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
//    bp::def("rhs_columnwise", blk_1m::rhs_columnwise); 
  }

  // blk_2m stuff
  {
    std::string nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".blk_2m");
    bp::object nested_module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
    bp::scope().attr("blk_2m") = nested_module;
    bp::scope parent = nested_module;
    bp::class_<libcloudphxx::blk_2m::opts_t<real_t>>("opts_t");
//    bp::def("rhs_cellwise", blk_2m::rhs_cellwise);
//    bp::def("rhs_columnwise", blk_2m::rhs_columnwise);
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
