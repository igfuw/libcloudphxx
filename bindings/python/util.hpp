// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

#pragma once

#if defined(BZ_THREADSAFE)
#  error please unset BZ_THREADSAFE
#endif
#include <blitz/tv2fastiter.h> // otherwise Clang fails in debug mode
#include <blitz/array.h>

#include <boost/python.hpp>

#include <libcloudph++/lgrngn/arrinfo.hpp>

namespace libcloudphxx
{
  namespace python
  {
    namespace bp = boost::python;
    using py_ptr_t = long; // TODO: acquire it using some decltype()

    void sanity_checks(const bp::numeric::array &arg)
    {
      // assuring double precision
      if (std::string(bp::extract<std::string>(arg.attr("dtype").attr("name"))) != "float64")
	throw std::runtime_error("dtype=float64 required for all passed arrays");

      // assuring contiguous layout
      if (!bp::extract<bool>(arg.attr("flags").attr("c_contiguous")))
	throw std::runtime_error("contiguous memory layout required");
    }

    template <class arr_t>
    arr_t np2bz(const bp::numeric::array &arg)
    {
      sanity_checks(arg);

      // wrapping the data into a Blitz++ array to get STL-container-like functionality
      return arr_t(
	// pointer to the data
	reinterpret_cast<typename arr_t::T_numtype*>(
	  (py_ptr_t)bp::extract<py_ptr_t>(arg.attr("ctypes").attr("data")) 
	), 
	// length of the array (regardless of the original dimensionality, we do 1D)
	blitz::shape(bp::extract<long>(arg.attr("size"))), 
	// ensure Blitz++ does not try to free the memory when done
	blitz::neverDeleteData
      );
    }

    template <class real_t>
    lgrngn::arrinfo_t<real_t> np2ai(
      const bp::numeric::array &arg,
      const std::array<int, 3> &sz
    ) {
      sanity_checks(arg);

      // C++ array dimensionality
      const int n_dims = 
        (sz[0]  > 0 && sz[1]  > 0 && sz[2]  > 0) ? 3 :
        (sz[0]  > 0 && sz[1] == 0 && sz[2]  > 0) ? 2 :
        (sz[0] == 0 && sz[1] == 0 && sz[2]  > 0) ? 1 :
        (sz[0]  > 1 && sz[1] == 0 && sz[2] == 0) ? 1 :
        0;

      std::vector<ptrdiff_t> strides(std::max(1, n_dims));

      // getting original strides from NumPy
      for (int i = 0; i < bp::len(arg.attr("strides")); ++i)
	strides[i] = bp::extract<ptrdiff_t>(arg.attr("strides")[i]) / sizeof(real_t);

      // overriding 1-element strides for non-z dimensions 
      // to make it work with single-column profiles
      switch (n_dims) 
      {
        case 3: // 3D array in C++
          switch (bp::len(arg.attr("shape")))
          {
            case 3: // 3D array in Python
              break; 
            case 1: // 1D array in Python
              strides[2] = strides[0];
              strides[0] = strides[1] = 0;
              break;
            default: // else
              throw std::runtime_error("incompatible array size: 3D set-up accepts either 1D profiles or 3D arrays");
          }
          break;
        case 2: // 2D arrays in C++
          switch (bp::len(arg.attr("shape"))) // array dimensionality in Python
          {
            case 2: // 2D array in Python
              break;
            case 1: // 1D array in Python
              strides[1] = strides[0];
              strides[0] = 0;
              break;
            default: 
              throw std::runtime_error("incompatible array size: 2D set-up accepts either 1D profiles or 2D arrays");
          }
          break;
        case 1: // 1D arrays in C++
          if (bp::len(arg.attr("shape")) != 1)
            throw std::runtime_error("incompatible array size: 1D set-up accepts only 1D arrays");
          break;
        case 0: // parcel set-up
          if (bp::len(arg.attr("shape")) != 1)
            throw std::runtime_error("incompatible array size: 0D parcel set-up accepts only 1D arrays");
          if (bp::extract<int>(arg.attr("shape")[0]) != 1) 
	    throw std::runtime_error("incompatible array size: 0D parcel set-up accepts single-element arrays only");
          break;
        default: 
          assert(false);
      }

      // checking profile length 
/* TODO: 1D horizontal slab not taken into account - probably this check will not be possible here
      if (strides[0] == 0 || n_dims == 1)
      {
	if (bp::extract<int>(arg.attr("shape")[0]) != sz[2])
	  throw std::runtime_error("incompatible array size: expecting nz-element profile");
      }
*/

      // getting data pointer from NumPy and returning
      return lgrngn::arrinfo_t<real_t>(
	reinterpret_cast<real_t*>(
	  (py_ptr_t)bp::extract<py_ptr_t>(arg.attr("ctypes").attr("data"))
	),
	strides
      );
    }
  };
};
