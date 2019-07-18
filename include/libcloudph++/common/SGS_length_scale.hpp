#pragma once

namespace libcloudphxx
{
  namespace common
  {
    namespace SGS_length_scale
    {
      // lambda = dz as in SAM and UWLCM
      template <typename real_t>
      quantity<si::length, real_t> vertical(
        quantity<si::length, real_t> dx
      )
      { return quantity<si::length, real_t>(dx);}

      template <typename real_t>
      quantity<si::length, real_t> vertical(
        quantity<si::length, real_t> dx,
        quantity<si::length, real_t> dz
      )
      { return quantity<si::length, real_t>(dz);}

      template <typename real_t>
      quantity<si::length, real_t> vertical(
        quantity<si::length, real_t> dx,
        quantity<si::length, real_t> dy,
        quantity<si::length, real_t> dz
      )
      { return quantity<si::length, real_t>(dz);}

      // lambda = (dx*dy*dz)^(1/3)
      template <typename real_t>
      quantity<si::length, real_t> geometric_mean(
        quantity<si::length, real_t> dx
      )
      { return quantity<si::length, real_t>(dx);}

      template <typename real_t>
      quantity<si::length, real_t> geometric_mean(
        quantity<si::length, real_t> dx,
        quantity<si::length, real_t> dz
      )
      { return quantity<si::length, real_t>(sqrt(dx*dz));}

      template <typename real_t>
      quantity<si::length, real_t> geometric_mean(
        quantity<si::length, real_t> dx,
        quantity<si::length, real_t> dy,
        quantity<si::length, real_t> dz
      )
      { return quantity<si::length, real_t>(pow((dx/si::metres)*(dy/si::metres)*(dz/si::metres), real_t(1./3.)) * si::metres);}

      // lambda = (dx+dy+dz) / 3
      template <typename real_t>
      quantity<si::length, real_t> arithmetic_mean(
        quantity<si::length, real_t> dx
      )
      { return quantity<si::length, real_t>(dx);}

      template <typename real_t>
      quantity<si::length, real_t> arithmetic_mean(
        quantity<si::length, real_t> dx,
        quantity<si::length, real_t> dz
      )
      { return quantity<si::length, real_t>((dx+dz) / real_t(2));}

      template <typename real_t>
      quantity<si::length, real_t> arithmetic_mean(
        quantity<si::length, real_t> dx,
        quantity<si::length, real_t> dy,
        quantity<si::length, real_t> dz
      )
      { return quantity<si::length, real_t>((dx + dy + dz) / real_t(3));}
    };
  };
};
