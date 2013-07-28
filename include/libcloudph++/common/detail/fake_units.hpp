/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief an "empty" Boost.units API for compilers that cannot digest Boost.units
  */

#pragma once
#include <cmath>

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
	  inline quantity(real_t value) : value(value) {}
	  inline quantity() {}

          // operators
	  inline quantity operator*(const quantity &q) const { return quantity(value * q.value); };
	  inline quantity operator/(const quantity &q) const { return quantity(value / q.value); };
	  inline quantity operator-(const quantity &q) const { return quantity(value - q.value); };
	  inline quantity operator+(const quantity &q) const { return quantity(value + q.value); };

          // cast to real_t
          inline operator real_t() const { return value; }
	};

	template <typename qntt_t, typename real_t>
        inline quantity<qntt_t, real_t> pow(
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
	  inline quantity<qntt_t, real_t> operator*(const real_t &a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

	  template <typename real_t>
	  inline quantity<qntt_t, real_t> operator/(const real_t &a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

	  template <typename real_t>
	  inline quantity<qntt_t, real_t> operator*(quantity<qntt_t, real_t> a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

	  template <typename real_t>
	  inline quantity<qntt_t, real_t> operator/(quantity<qntt_t, real_t> a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

	  template <typename real_t>
	  inline quantity<qntt_t, real_t> operator*(const unit_t &, quantity<qntt_t, real_t> a) { return quantity<qntt_t, real_t>(a); }

	  //template <typename real_t>
	  //inline quantity<qntt_t, real_t> operator/(const unit_t &, quantity<qntt_t, real_t> a) { return pow(a, -1); }
	}

	namespace si
	{
	  typedef detail::qntt_t pressure;
	  typedef detail::qntt_t velocity;
	  typedef detail::qntt_t energy;
	  typedef detail::qntt_t mass_density;
	  typedef detail::qntt_t mass;
	  typedef detail::qntt_t temperature;
	  typedef detail::qntt_t dimensionless;
	  typedef detail::qntt_t amount;

	  static const detail::unit_t 
	    cubic_metre, cubic_metres,
	    joule, joules,
	    kelvin, kelvins,
	    kilogram, kilograms,
	    metre, metres, 
	    mole, moles,
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
