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
          return quantity<qntt_t, real_t>(std::pow(a.value, b.value));
        }

	// ...
	namespace detail
	{
	  struct qntt_t {};
     
	  struct unit_t 
	  {
	    unit_t() {}
	  };

	  template <typename real_t>
	  BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator*(const real_t &a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

	  template <typename real_t>
	  BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator/(const real_t &a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

	  template <typename real_t>
	  BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator*(quantity<qntt_t, real_t> a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

	  template <typename real_t>
	  BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator/(quantity<qntt_t, real_t> a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

	  template <typename real_t>
	  BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator*(const unit_t &, quantity<qntt_t, real_t> a) { return quantity<qntt_t, real_t>(a); }

	  template <typename real_t>
	  BOOST_GPU_ENABLED inline quantity<qntt_t, real_t> operator+(const int &c, quantity<qntt_t, real_t> a) { return a+c; }
	}

	namespace si
	{
	  typedef detail::qntt_t amount;
	  typedef detail::qntt_t dimensionless;
	  typedef detail::qntt_t energy;
	  typedef detail::qntt_t force;
	  typedef detail::qntt_t length;
	  typedef detail::qntt_t mass;
	  typedef detail::qntt_t mass_density;
	  typedef detail::qntt_t pressure;
	  typedef detail::qntt_t temperature;
	  typedef detail::qntt_t velocity;
	  typedef detail::qntt_t volume;

	  static const detail::unit_t 
	    cubic_metre, cubic_metres,
	    joule, joules,
	    kelvin, kelvins,
	    kilogram, kilograms,
	    metre, metres, 
	    mole, moles,
            newton, newtons,
	    pascal, pascals;
        };

        template <typename, typename>    
        struct divide_typeof_helper { typedef detail::qntt_t type; };

        template <typename, typename>    
        struct multiply_typeof_helper { typedef detail::qntt_t type; };
      };
    };
  };
};
