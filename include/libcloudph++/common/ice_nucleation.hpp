#pragma once

#include "units.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace ice_nucleation
    {
      enum class INP_t {mineral, soot}; // types of ice nucleating particles

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<divide_typeof_helper<si::dimensionless, si::area>::type, real_t> n_s(
        const quantity<si::temperature, real_t> T,        // ambient temperature
        const INP_t& INP_type                             // type of ice nucleating particle
      ) {

        if (INP_type == INP_t::mineral) // valid between 237 - 261 K
        {
          return std::exp(-real_t(0.517) * (T/si::kelvins - real_t(273.15)) + real_t(8.934)) / si::square_meters;
        }
        if (INP_type == INP_t::soot) // valid between 239 - 255 K
        {
          return real_t(7.463)*std::exp(-real_t(0.0101) * std::pow(T/si::kelvins - real_t(273.15) , 2) - real_t(0.8525) * (T/si::kelvins - real_t(273.15)) + real_t(0.7667)) / si::square_meters;
        }
      }

      using dn_dT_dimension = boost::units::divide_typeof_helper<
        boost::units::divide_typeof_helper<boost::units::si::dimensionless, boost::units::si::area>::type,
        boost::units::si::temperature
      >::type; // m^-2 / K

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<dn_dT_dimension, real_t> dn_dT(
        const quantity<si::temperature, real_t> T,        // ambient temperature
        const INP_t& INP_type                             // type of ice nucleating particle
      ) {

        if (INP_type == INP_t::mineral) // valid between 237 - 261 K
        {
          return -real_t(0.517) * n_s(T, INP_type) / si::kelvins;
        }
        if (INP_type == INP_t::soot) // valid between 239 - 255 K
        {
          return (-real_t(0.0202) * T/si::kelvins - real_t(0.8525)) * n_s(T, INP_type) / si::kelvins;
        }
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<divide_typeof_helper<si::dimensionless, si::temperature>::type, real_t> p(
        const quantity<si::temperature, real_t> T,      // temperature
        const INP_t& INP_type,                          // type of ice nucleating particle
        const quantity<si::length, real_t> rd_insol     // radius of ice nucleating (insoluble) particle
      ) {
        const quantity<si::area, real_t> A = real_t(4) * pi<real_t>() * rd_insol * rd_insol; // surface area of the insoluble particle
        return -A * dn_dT(T, INP_type) * std::exp(-A * n_s(T, INP_type));
      }

    };
  };
};
