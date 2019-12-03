/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief an "empty" Boost.units API for compilers that cannot digest Boost.units
  */

#pragma once
#include <cmath> // std::pow
#include <boost/config.hpp> // BOOST_GPU_ENABLED

// TODO: CUDA's _forceinline_?
// TODO: CUDA's _device_

namespace libcloudphxx
{
  namespace common
  {
    namespace detail
    {
      namespace fake_units
      {
      template <typename, typename real_t>
      struct quantity
      {
        real_t value;

          // ctors
        BOOST_GPU_ENABLED inline quantity(real_t value) : value(value) {}
        BOOST_GPU_ENABLED inline quantity() {}

          // operators
        BOOST_GPU_ENABLED inline quantity operator*(const quantity &q) const { return quantity(value * q.value); };
        BOOST_GPU_ENABLED inline quantity operator/(const quantity &q) const { return quantity(value / q.value); };
        BOOST_GPU_ENABLED inline quantity operator-(const quantity &q) const { return quantity(value - q.value); };
        BOOST_GPU_ENABLED inline quantity operator+(const quantity &q) const { return quantity(value + q.value); };

        BOOST_GPU_ENABLED inline quantity operator+(const int         &q) const { return quantity(value + q); };

          // cast to real_t
          BOOST_GPU_ENABLED inline operator real_t() const { return value; }
      };

        // pow function
      template <typename qntt_t, typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> pow(
          const quantity<qntt_t, real_t> &a,
          const quantity<qntt_t, real_t> &b
        )
        {
          return quantity<qntt_t, real_t>(pow(a.value, b.value));
        }

      // ...
      namespace detail
      {
        struct qntt_t {};
     
        struct unit_t {
            unit_t() {}
          };

          // real_t vs. unit
        template <typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator*(const real_t &a, const unit_t &) { return quantity<qntt_t, real_t>(a); }
        template <typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator/(const real_t &a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

          // quantity vs. unit
        template <typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator*(quantity<qntt_t, real_t> a, const unit_t &) { return quantity<qntt_t, real_t>(a); }
        template <typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator/(quantity<qntt_t, real_t> a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

          // unit vs. quantity
        template <typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator*(const unit_t &, quantity<qntt_t, real_t> a) { return quantity<qntt_t, real_t>(a); }

          // num_t vs. quantity
        template <typename num_t, typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator+(const num_t &c, quantity<qntt_t, real_t> a) { return quantity<qntt_t, real_t>(c)+a; }
        template <typename num_t, typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator-(const num_t &c, quantity<qntt_t, real_t> a) { return quantity<qntt_t, real_t>(c)-a; }
        template <typename num_t, typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator/(const num_t &c, quantity<qntt_t, real_t> a) { return quantity<qntt_t, real_t>(c)/a; }
        template <typename num_t, typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator*(const num_t &c, quantity<qntt_t, real_t> a) { return quantity<qntt_t, real_t>(c)*a; }

          // quantity vs. num_t
        template <typename num_t, typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator+(quantity<qntt_t, real_t> a, const num_t &c) { return a+quantity<qntt_t, real_t>(c); }
        template <typename num_t, typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator-(quantity<qntt_t, real_t> a, const num_t &c) { return a-quantity<qntt_t, real_t>(c); }
        template <typename num_t, typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator/(quantity<qntt_t, real_t> a, const num_t &c) { return a/quantity<qntt_t, real_t>(c); }
        template <typename num_t, typename real_t>
        BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator*(quantity<qntt_t, real_t> a, const num_t &c) { return a*quantity<qntt_t, real_t>(c); }
      }

      namespace si
      {
        typedef detail::qntt_t acceleration;
        typedef detail::qntt_t amount;
        typedef detail::qntt_t area;
        typedef detail::qntt_t dimensionless;
        typedef detail::qntt_t dynamic_viscosity;
        typedef detail::qntt_t energy;
        typedef detail::qntt_t force;
        typedef detail::qntt_t length;
        typedef detail::qntt_t mass;
        typedef detail::qntt_t mass_density;
        typedef detail::qntt_t pressure;
        typedef detail::qntt_t temperature;
        typedef detail::qntt_t time;
        typedef detail::qntt_t velocity;
        typedef detail::qntt_t volume;

// if needed since this file is included in cpp test fake_units
#if defined(__NVCC__)
          __device__
#endif
        static const detail::unit_t 
          cubic_metre, cubic_metres, cubic_meter, cubic_meters,
          joule, joules,
          kelvin, kelvins,
          kilogram, kilograms,
          metre, metres, meter, meters, 
          metres_per_second, meters_per_second, 
          metres_per_second_squared, meters_per_second_squared, 
          mole, moles,
            newton, newtons,
          pascal, pascals,
            second, seconds,
            square_metre, square_metres, square_meter, square_meters;
        };

        template <typename, typename>    
        struct divide_typeof_helper { typedef detail::qntt_t type; };

        template <typename, typename>    
        struct multiply_typeof_helper { typedef detail::qntt_t type; };

        template <typename, typename>
        struct power_typeof_helper { typedef detail::qntt_t type; };

        template <int>
        struct static_rational {};
      };
    };
  };
};
