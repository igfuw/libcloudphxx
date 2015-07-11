//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// adapted for use on GPU by Maciej Waruszewski

#pragma once

#include <boost/math/tools/config.hpp>

#include <cfloat>
#include <cassert>

namespace libcloudphxx {
namespace common {
namespace detail {
namespace toms748_detail
{
  namespace // anonymous namespace to prevent multiple definition error at link time
  {
    template<class T> T max_value();
    template<> BOOST_GPU_ENABLED float max_value<float>() {return FLT_MAX;}
    template<> BOOST_GPU_ENABLED double max_value<double>() {return DBL_MAX;}

    template<class T> T min_value();
    template<> BOOST_GPU_ENABLED float min_value<float>() {return FLT_MIN;}
    template<> BOOST_GPU_ENABLED double min_value<double>() {return DBL_MIN;}

    template<class T> T epsilon();
    template<> BOOST_GPU_ENABLED float epsilon<float>() {return FLT_EPSILON;}
    template<> BOOST_GPU_ENABLED double epsilon<double>() {return DBL_EPSILON;}
  }

  template <class T>
  class eps_tolerance
  {
  public:
     BOOST_GPU_ENABLED
     eps_tolerance(unsigned bits)
     {
#if !defined(__NVCC__)
        using std::max;
#endif
        eps = max(T(ldexp(1.0F, 1-bits)), T(4 * epsilon<T>()));
     }
     BOOST_GPU_ENABLED
     bool operator()(const T& a, const T& b)
     {
#if !defined(__NVCC__)
        using std::min;
#endif
        return fabs(a - b) <= (eps * min(fabs(a), fabs(b)));
     }
  private:
     T eps;
  };

  template <class F, class T>
  BOOST_GPU_ENABLED
  void bracket(F f, T& a, T& b, T c, T& fa, T& fb, T& d, T& fd)
  {
     //
     // Given a point c inside the existing enclosing interval
     // [a, b] sets a = c if f(c) == 0, otherwise finds the new 
     // enclosing interval: either [a, c] or [c, b] and sets
     // d and fd to the point that has just been removed from
     // the interval.  In other words d is the third best guess
     // to the root.
     //
     T tol = epsilon<T>() * 2;
     //
     // If the interval [a,b] is very small, or if c is too close 
     // to one end of the interval then we need to adjust the
     // location of c accordingly:
     //
     if((b - a) < 2 * tol * a)
     {
        c = a + (b - a) / 2;
     }
     else if(c <= a + fabs(a) * tol)
     {
        c = a + fabs(a) * tol;
     }
     else if(c >= b - fabs(b) * tol)
     {
        c = b - fabs(a) * tol;
     }
     //
     // OK, lets invoke f(c):
     //
     T fc = f(c);
     //
     // if we have a zero then we have an exact solution to the root:
     //
     if(fc == 0)
     {
        a = c;
        fa = 0;
        d = 0;
        fd = 0;
        return;
     }
     //
     // Non-zero fc, update the interval:
     //
     if(copysign(T(1), fa * fc) < 0)
     {
        d = b;
        fd = fb;
        b = c;
        fb = fc;
     }
     else
     {
        d = a;
        fd = fa;
        a = c;
        fa= fc;
     }
  }

  template <class T>
  BOOST_GPU_ENABLED
  inline T safe_div(T num, T denom, T r)
  {
     //
     // return num / denom without overflow,
     // return r if overflow would occur.
     //

     if(fabs(denom) < 1)
     {
        if(fabs(denom * max_value<T>()) <= fabs(num))
           return r;
     }
     return num / denom;
  }

  template <class T>
  BOOST_GPU_ENABLED
  inline T secant_interpolate(const T& a, const T& b, const T& fa, const T& fb)
  {
     //
     // Performs standard secant interpolation of [a,b] given
     // function evaluations f(a) and f(b).  Performs a bisection
     // if secant interpolation would leave us very close to either
     // a or b.  Rationale: we only call this function when at least
     // one other form of interpolation has already failed, so we know
     // that the function is unlikely to be smooth with a root very
     // close to a or b.
     //

     T tol = epsilon<T>() * 5;
     T c = a - (fa / (fb - fa)) * (b - a);
     if((c <= a + fabs(a) * tol) || (c >= b - fabs(b) * tol))
        return (a + b) / 2;
     return c;
  }

  template <class T>
  BOOST_GPU_ENABLED
  T quadratic_interpolate(const T& a, const T& b, T const& d,
                             const T& fa, const T& fb, T const& fd, 
                             unsigned count)
  {
     //
     // Performs quadratic interpolation to determine the next point,
     // takes count Newton steps to find the location of the
     // quadratic polynomial.
     //
     // Point d must lie outside of the interval [a,b], it is the third
     // best approximation to the root, after a and b.
     //
     // Note: this does not guarentee to find a root
     // inside [a, b], so we fall back to a secant step should
     // the result be out of range.
     //
     // Start by obtaining the coefficients of the quadratic polynomial:
     //
     T B = safe_div(T(fb - fa), T(b - a), max_value<T>());
     T A = safe_div(T(fd - fb), T(d - b), max_value<T>());
     A = safe_div(T(A - B), T(d - a), T(0));

     if(A == 0)
     {
        // failure to determine coefficients, try a secant step:
        return secant_interpolate(a, b, fa, fb);
     }
     //
     // Determine the starting point of the Newton steps:
     //
     T c;
     if(copysign(T(1), A * fa) > 0)
     {
        c = a;
     }
     else
     {
        c = b;
     }
     //
     // Take the Newton steps:
     //
     for(unsigned i = 1; i <= count; ++i)
     {
        //c -= safe_div(B * c, (B + A * (2 * c - a - b)), 1 + c - a);
        c -= safe_div(T(fa+(B+A*(c-b))*(c-a)), T(B + A * (2 * c - a - b)), T(1 + c - a));
     }
     if((c <= a) || (c >= b))
     {
        // Oops, failure, try a secant step:
        c = secant_interpolate(a, b, fa, fb);
     }
     return c;
  }

  template <class T>
  BOOST_GPU_ENABLED
  T cubic_interpolate(const T& a, const T& b, const T& d, 
                      const T& e, const T& fa, const T& fb, 
                      const T& fd, const T& fe)
  {
     //
     // Uses inverse cubic interpolation of f(x) at points 
     // [a,b,d,e] to obtain an approximate root of f(x).
     // Points d and e lie outside the interval [a,b]
     // and are the third and forth best approximations
     // to the root that we have found so far.
     //
     // Note: this does not guarentee to find a root
     // inside [a, b], so we fall back to quadratic
     // interpolation in case of an erroneous result.
     //


     T q11 = (d - e) * fd / (fe - fd);
     T q21 = (b - d) * fb / (fd - fb);
     T q31 = (a - b) * fa / (fb - fa);
     T d21 = (b - d) * fd / (fd - fb);
     T d31 = (a - b) * fb / (fb - fa);

     T q22 = (d21 - q11) * fb / (fe - fb);
     T q32 = (d31 - q21) * fa / (fd - fa);
     T d32 = (d31 - q21) * fd / (fd - fa);
     T q33 = (d32 - q22) * fa / (fe - fa);
     T c = q31 + q32 + q33 + a;

     if((c <= a) || (c >= b))
     {
        // Out of bounds step, fall back to quadratic interpolation:
        c = quadratic_interpolate(a, b, d, fa, fb, fd, 3);
     }

     return c;
  }
} // namespace toms748_detail

template <class T>
class eps_tolerance
{
  public:

  BOOST_GPU_ENABLED
  eps_tolerance(unsigned bits)
  {
#if !defined(__NVCC__)
    using std::max;
#endif
    eps = max(T(ldexp(1.0F, 1-bits)), T(4 * toms748_detail::epsilon<T>()));
  }

  BOOST_GPU_ENABLED
  bool operator()(const T& a, const T& b)
  {
#if !defined(__NVCC__)
    using std::min;
#endif
    return fabs(a - b) <= (eps * min(fabs(a), fabs(b)));
  }

  private:
  T eps;
};


template <class F, class T, class Tol>
BOOST_GPU_ENABLED
T toms748_solve(F f, const T& ax, const T& bx, const T& fax, const T& fbx, Tol tol, uintmax_t &max_iter)
{
   uintmax_t count = max_iter;
   T a, b, fa, fb, c, u, fu, a0, b0, d, fd, e, fe;
   const T mu = 0.5f;

   a = ax;
   b = bx;
   
   assert(a < b);

   fa = fax;
   fb = fbx;

   if(tol(a, b) || (fa == 0) || (fb == 0))
   {
      max_iter = 0;
      if(fa == 0)
         b = a;
      else if(fb == 0)
         a = b;
      return (a + b)/2;
   }

   assert(copysign(T(1), fa * fb) < 0);

   // dummy value for fd, e and fe:
   fe = e = fd = 1e5F;

   if(fa != 0)
   {
      //
      // On the first step we take a secant step:
      //
      c = toms748_detail::secant_interpolate(a, b, fa, fb);
      toms748_detail::bracket(f, a, b, c, fa, fb, d, fd);
      --count;

      if(count && (fa != 0) && !tol(a, b))
      {
         //
         // On the second step we take a quadratic interpolation:
         //
         c = toms748_detail::quadratic_interpolate(a, b, d, fa, fb, fd, 2);
         e = d;
         fe = fd;
         toms748_detail::bracket(f, a, b, c, fa, fb, d, fd);
         --count;
      }
   }

   while(count && (fa != 0) && !tol(a, b))
   {
      // save our brackets:
      a0 = a;
      b0 = b;
      //
      // Starting with the third step taken
      // we can use either quadratic or cubic interpolation.
      // Cubic interpolation requires that all four function values
      // fa, fb, fd, and fe are distinct, should that not be the case
      // then variable prof will get set to true, and we'll end up
      // taking a quadratic step instead.
      //
      T min_diff = toms748_detail::min_value<T>() * 32;
      bool prof = (fabs(fa - fb) < min_diff) || (fabs(fa - fd) < min_diff) || (fabs(fa - fe) < min_diff) || (fabs(fb - fd) < min_diff) || (fabs(fb - fe) < min_diff) || (fabs(fd - fe) < min_diff);
      if(prof)
      {
         c = toms748_detail::quadratic_interpolate(a, b, d, fa, fb, fd, 2);
      }
      else
      {
         c = toms748_detail::cubic_interpolate(a, b, d, e, fa, fb, fd, fe);
      }
      //
      // re-bracket, and check for termination:
      //
      e = d;
      fe = fd;
      toms748_detail::bracket(f, a, b, c, fa, fb, d, fd);
      if((0 == --count) || (fa == 0) || tol(a, b))
         break;
      //
      // Now another interpolated step:
      //
      prof = (fabs(fa - fb) < min_diff) || (fabs(fa - fd) < min_diff) || (fabs(fa - fe) < min_diff) || (fabs(fb - fd) < min_diff) || (fabs(fb - fe) < min_diff) || (fabs(fd - fe) < min_diff);
      if(prof)
      {
         c = toms748_detail::quadratic_interpolate(a, b, d, fa, fb, fd, 3);
      }
      else
      {
         c = toms748_detail::cubic_interpolate(a, b, d, e, fa, fb, fd, fe);
      }
      //
      // Bracket again, and check termination condition, update e:
      //
      toms748_detail::bracket(f, a, b, c, fa, fb, d, fd);
      if((0 == --count) || (fa == 0) || tol(a, b))
         break;
      //
      // Now we take a double-length secant step:
      //
      if(fabs(fa) < fabs(fb))
      {
         u = a;
         fu = fa;
      }
      else
      {
         u = b;
         fu = fb;
      }
      c = u - 2 * (fu / (fb - fa)) * (b - a);
      if(fabs(c - u) > (b - a) / 2)
      {
         c = a + (b - a) / 2;
      }
      //
      // Bracket again, and check termination condition:
      //
      e = d;
      fe = fd;
      toms748_detail::bracket(f, a, b, c, fa, fb, d, fd);
      if((0 == --count) || (fa == 0) || tol(a, b))
         break;
      //
      // And finally... check to see if an additional bisection step is 
      // to be taken, we do this if we're not converging fast enough:
      //
      if((b - a) < mu * (b0 - a0))
         continue;
      //
      // bracket again on a bisection:
      //
      e = d;
      fe = fd;
      toms748_detail::bracket(f, a, b, T(a + (b - a) / 2), fa, fb, d, fd);
      --count;
   } // while loop

   max_iter -= count;
   if(fa == 0)
   {
      b = a;
   }
   else if(fb == 0)
   {
      a = b;
   }
   return (a + b)/2;
}

template <class F, class T, class Tol>
BOOST_GPU_ENABLED
inline T toms748_solve(F f, const T& ax, const T& bx, Tol tol, uintmax_t &max_iter)
{
   max_iter -= 2;
   T r = toms748_solve(f, ax, bx, f(ax), f(bx), tol, max_iter);
   max_iter += 2;
   return r;
}

// versions with hopefully-sensible deafult tolerance and number of iterations

template <class F, class T>
BOOST_GPU_ENABLED
T toms748_solve(F f, const T& ax, const T& bx, const T& fax, const T& fbx)
{
  const uintmax_t n_iter = 20;
  uintmax_t max_iter = n_iter;
  T root = toms748_solve(f, ax, bx, fax, fbx, 
    common::detail::eps_tolerance<T>(sizeof(T) * 8 / 4),
    max_iter
  );
  assert(max_iter != n_iter);
  return root;
}

template <class F, class T>
BOOST_GPU_ENABLED
T toms748_solve(F f, const T& ax, const T& bx)
{
  return toms748_solve(f, ax, bx, f(ax), f(bx));
}

} // namespace detail
} // namespace common
} // namespace libcloudphxx
