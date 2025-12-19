// Python binding for libcloudph++
//
// author[s]: Sylwester Arabas, Dorota Jarecka
// licensing: GPU GPL v3
// copyright: University of Warsaw

#include "error.hpp"

#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/kappa_koehler.hpp>
#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/henry.hpp>
#include <libcloudph++/common/dissoc.hpp>
#include <libcloudph++/common/tetens.hpp>

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
      real_t exner(const real_t &p)
      {
	return cmn::theta_std::exner(p * si::pascals);
      }

      template <typename real_t>
      real_t p_v(const real_t &p, const real_t &r)
      {
	return cmn::moist_air::p_v(p * si::pascals, r * si::dimensionless()) / si::pascals;
      }

      template <typename real_t>
      real_t p_vs(const real_t &T)
      {
	return cmn::const_cp::p_vs(T * si::kelvins) / si::pascals;
      }

      template <typename real_t>
      real_t p_vsi(const real_t &T)
      {
        return cmn::const_cp::p_vsi(T * si::kelvins) / si::pascals;
      }

      template <typename real_t>
      real_t p_vs_tet(const real_t &T)
      {
	return cmn::tetens::p_vs(T * si::kelvins) / si::pascals;
      }
      
      template <typename real_t>
      real_t r_vs(const real_t &T, const real_t &p)
      {
	return cmn::const_cp::r_vs(T * si::kelvins, p * si::pascals);
      }
      
      template <typename real_t>
      real_t r_vsi(const real_t &T, const real_t &p)
      {
	return cmn::const_cp::r_vsi(T * si::kelvins, p * si::pascals);
      }
	    
      template <typename real_t>
      real_t l_v(const real_t &T)
      {
	return cmn::const_cp::l_v(T * si::kelvins) / si::joules * si::kilograms;
      }

      template <typename real_t>
      real_t l_s(const real_t &T)
      {
        return cmn::const_cp::l_s(T * si::kelvins) / si::joules * si::kilograms;
      }

      template <typename real_t>
      real_t l_f(const real_t &T)
      {
        return cmn::const_cp::l_f(T * si::kelvins) / si::joules * si::kilograms;
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
      real_t visc(const real_t &T)
      {
	return cmn::vterm::visc(T * si::kelvins) / si::pascals / si::seconds;
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

      template <typename real_t>
      real_t p_hydro(
        const real_t &z,
        const real_t &th_0,
        const real_t &r_0,
        const real_t &z_0,
        const real_t &p_0
      )
      {
        return cmn::hydrostatic::p(
          z * si::metres, 
          th_0 * si::kelvins, 
          r_0 * si::dimensionless(), 
          z_0 * si::metres, 
          p_0 * si::pascals
        ) / si::pascals;
      }

      template <typename real_t>
      real_t rhod(
        const real_t &p,
        const real_t &th_std,
        const real_t &r_v
      )
      {
        return cmn::theta_std::rhod(
          p * si::pascals,
          th_std * si::kelvins,
          r_v * si::dimensionless()
        ) / si::kilograms * si::cubic_metres;
      }
    };
  };
};
