// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

#include "error.hpp"

#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/kappa_koehler.hpp>

namespace libcloudphxx
{
  namespace python
  {
    namespace cmn = libcloudphxx::common;

    namespace common
    {
      template <typename real_t>
      real_t th_dry2std(
	const real_t &th_dry,
	const real_t &r
      )
      {
	return cmn::theta_dry::dry2std(th_dry * si::kelvins, r * si::dimensionless()) / si::kelvins;
      }

      template <typename real_t>
      real_t th_std2dry(
	const real_t &th_std,
	const real_t &r
      )
      {
	return cmn::theta_dry::std2dry(th_std * si::kelvins, r * si::dimensionless()) / si::kelvins;
      }

      template <typename real_t>
      real_t p_vs(const real_t &T)
      {
	return cmn::const_cp::p_vs(T * si::kelvins) / si::pascals;
      }

      template <typename real_t>
      real_t T(const real_t &th, const real_t &rhod)
      {
	return cmn::theta_dry::T(th * si::kelvins, rhod  * si::kilograms / si::cubic_metres) / si::kelvins;
      }

      template <typename real_t>
      real_t p(const real_t &rhod, const real_t &r, const real_t &T)
      {
	return cmn::theta_dry::p(rhod  * si::kilograms / si::cubic_metres, r * si::kilograms / si::kilograms, T * si::kelvins) / si::pascals;
      }

      template <typename real_t>
      real_t rw3_cr(const real_t &rd3, const real_t &kappa, const real_t &T)
      {
        return cmn::kappa_koehler::rw3_cr(
          rd3 * si::cubic_metres,
          kappa * si::dimensionless(),
          T * si::kelvins
        ) / si::cubic_metres;
      }

      template <typename real_t>
      real_t S_cr(const real_t &rd3, const real_t &kappa, const real_t &T)
      {
        return cmn::kappa_koehler::S_cr(
          rd3 * si::cubic_metres,
          kappa * si::dimensionless(),
          T * si::kelvins
        );
      }
    };
  };
};
