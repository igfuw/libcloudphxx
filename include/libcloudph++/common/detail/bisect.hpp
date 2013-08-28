#pragma once
// adapten from Boost's bisection

#include <boost/config.hpp>

#if !defined(__NVCC__)
#  include <cmath>
#  include <cassert>
#endif

#include <iostream>

namespace libcloudphxx
{
  namespace common
  {
    namespace detail
    {

      template <class f_t, typename real_t>
      BOOST_GPU_ENABLED
      inline real_t bisect(
        const f_t f, 
        real_t min, real_t max, 
        const real_t &tol,
        real_t fmin, real_t fmax
      )
      {
        assert(min < max); // assumes in order

        if (fmin == 0) return min;
        if (fmax == 0) return max;

#if !defined(__NVCC__)
	using std::abs;
#endif
        if (fmin * fmax >= 0) return (min + max) / 2; 

	real_t mid; 
	while (
          mid = (min + max) / 2, 
          abs(max - min) > tol
        )
	{
          if (mid == max || mid == min) break;

          real_t fmid = f(mid);

          if (fmid == 0) break;
          else if (fmid * fmin > 0) 
	    min = mid, fmin = fmid;
	  else 
	    max = mid, fmax = fmid;
	}
	return mid;
      }

      template <class f_t, typename real_t>
      BOOST_GPU_ENABLED
      inline real_t bisect(
        const f_t f, 
        const real_t &min, const real_t &max, 
        const real_t &tol,
        const real_t &fmin
      )
      {
        return bisect(f, min, max, tol, fmin, f(max));
      }

      template <class f_t, typename real_t>
      BOOST_GPU_ENABLED
      inline real_t bisect(
        const f_t f, 
        const real_t &min, const real_t &max, 
        const real_t &tol
      )
      {
        return bisect(f, min, max, tol, f(min), f(max));
      }
    };
  };
};
