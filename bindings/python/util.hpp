// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#if defined(BZ_THREADSAFE)
#  error please unset BZ_THREADSAFE
#endif
#include <blitz/array.h>

#include <boost/python.hpp>

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

inline lgr::arrinfo_t<real_t> np2ai(const bp::numeric::array &arg)
{
  sanity_checks(arg);

  const ptrdiff_t one = 1;
  
  return lgr::arrinfo_t<real_t>(
    reinterpret_cast<real_t*>(
      (py_ptr_t)bp::extract<py_ptr_t>(arg.attr("ctypes").attr("data"))
    ),
    &one // TODO: parcel assumption hardcoded
  );
}
