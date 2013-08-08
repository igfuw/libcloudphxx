#pragma once
// based on http://pegasus.rutgers.edu/~randall/473fall2004/hand02.pdf

#include <boost/config.hpp>

#if !defined(__NVCC__)
#  include <cmath>
#  include <cassert>
#endif

namespace libcloudphxx
{
  namespace common
  {
    namespace detail
    {
      template <class f_t, typename real_t>
      BOOST_GPU_ENABLED
      real_t bisect(f_t f, real_t a, real_t b, real_t tolerance)
      {
      #if !defined(__NVCC__)
	assert(f(a) * f(b) < 0); // must have different signs
	using std::abs;
	using std::sqrt;
      #endif

	real_t m;
	while (abs(b-a) > tolerance)
	{
	  m = (a+b) / 2;
	  real_t ae = abs(sqrt(2.0) - m); // actual error
	  real_t me = abs((b-a)/2); // max error
	  if (f(a) * f(m) > 0) 
	    a = m; 
	  else 
	    b = m;
	}
	return m;
      }
    };
  };
};
