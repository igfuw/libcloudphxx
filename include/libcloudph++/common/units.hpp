/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Common include directives and using statements for Boost.units
  */

#pragma once

#if !defined(NVCC)

#  include <boost/units/systems/si.hpp>
#  include <boost/units/cmath.hpp>
#  include <boost/units/io.hpp>
   namespace si = boost::units::si;
   using boost::units::quantity;
   using boost::units::one;
   using boost::units::pow;
   using boost::units::root;
   using boost::units::multiply_typeof_helper;
   using boost::units::divide_typeof_helper;
   using boost::units::power_typeof_helper;
   using boost::units::static_rational;

#  include <boost/math/constants/constants.hpp>
   using boost::math::constants::pi;

#else

   template <typename, typename real_t>
   struct quantity
   {
     real_t value;
     inline quantity(real_t value) : value(value) {}
     //quantity& operator=(const quantity &q) { value = q; };
   };

   // ...
   namespace fake
   {
     struct qntt_t {};

     struct unit_t 
     {
     };

     template <typename real_t>
     inline quantity<qntt_t, real_t> operator*(real_t &a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

     template <typename real_t>
     inline quantity<qntt_t, real_t> operator/(real_t &a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

/* TODO... does not work yet :(
     template <typename real_t>
     inline quantity<qntt_t, real_t> operator*(quantity<qntt_t, real_t> &a, const unit_t &) { return quantity<qntt_t, real_t>(a); }

     template <typename real_t>
     inline quantity<qntt_t, real_t> operator/(quantity<qntt_t, real_t> &a, const unit_t &) { return quantity<qntt_t, real_t>(a); }
*/
   }

   namespace si
   {
     typedef fake::qntt_t pressure;
     typedef fake::qntt_t velocity;
     typedef fake::qntt_t energy;
     typedef fake::qntt_t mass_density;
     typedef fake::qntt_t mass;
     typedef fake::qntt_t temperature;
     typedef fake::qntt_t dimensionless;
     typedef fake::qntt_t amount;

     static const fake::unit_t pascals;
     static const fake::unit_t metres;
     static const fake::unit_t cubic_metres;
     static const fake::unit_t kilograms;
     static const fake::unit_t kelvins;
   };

   template <typename, typename>    
   struct divide_typeof_helper { typedef fake::qntt_t type; };

   template <typename, typename>    
   struct multiply_typeof_helper { typedef fake::qntt_t type; };

#endif
