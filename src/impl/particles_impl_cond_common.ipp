// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/iterator/transform_iterator.h>
#include <libcloudph++/common/maxwell-mason.hpp>
#include <libcloudph++/common/kappa_koehler.hpp>
#include <libcloudph++/common/kelvin_term.hpp>
#include <libcloudph++/common/transition_regime.hpp>
#include <libcloudph++/common/ventil.hpp>
#include <libcloudph++/common/mean_free_path.hpp>
#include <libcloudph++/common/detail/toms748.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class real_t>
      struct rw3diff2drv 
      {
        real_t mlt;
        int n_dims;
        rw3diff2drv(const real_t &mlt, const int &n_dims):
          mlt(mlt), n_dims(n_dims) {}
 
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw3diff, const thrust::tuple<real_t, real_t, real_t> &tpl)
        {
          if(n_dims > 0)
            return mlt * rw3diff * thrust::get<1>(tpl) / thrust::get<0>(tpl) / thrust::get<2>(tpl);
          else // for parcel setup use 1/rhod instead of dv, dv will be updated in hskpng_Tpr in async
            return mlt * rw3diff * thrust::get<1>(tpl);
        }
      };

      template<class real_t>
      struct rw2torw3 : thrust::unary_function<const real_t&, real_t>
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2)
        {
#if !defined(__NVCC__)
          using std::sqrt;
#endif
          const real_t rw = sqrt(rw2);
          return rw2 * rw;
        }
      };
        
      template <typename real_t>
      struct advance_rw2_minfun
      {
        const quantity<si::area,              real_t> rw2_old;
        const quantity<si::time,              real_t> dt;
        const quantity<si::mass_density,      real_t> rhod;
        const quantity<si::dimensionless,     real_t> rv;
        const quantity<si::temperature,       real_t> T;
        const quantity<si::pressure,          real_t> p;
        const quantity<si::dimensionless,     real_t> RH;
        const quantity<si::dynamic_viscosity, real_t> eta;
        const quantity<si::volume,            real_t> rd3;
        const quantity<si::dimensionless,     real_t> kpa;
        const quantity<si::velocity,          real_t> vt;
        const quantity<si::dimensionless,     real_t> RH_max;
        const quantity<si::length,            real_t> lambda_D;
        const quantity<si::length,            real_t> lambda_K;

        // ctor
        BOOST_GPU_ENABLED
        advance_rw2_minfun(
          const real_t &dt,
          const real_t &rw2,
          const thrust::tuple<thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t>, real_t, real_t> &tpl,
          const real_t &RH_max
        ) : 
          dt(dt * si::seconds), 
          rw2_old(rw2 * si::square_metres),
          rhod(    thrust::get<0>(thrust::get<0>(tpl)) * si::kilograms / si::cubic_metres),
          rv(      thrust::get<1>(thrust::get<0>(tpl))),
          T(       thrust::get<2>(thrust::get<0>(tpl)) * si::kelvins),
          eta(     thrust::get<3>(thrust::get<0>(tpl)) * si::pascals * si::seconds),
          rd3(     thrust::get<4>(thrust::get<0>(tpl)) * si::cubic_metres),
          kpa(     thrust::get<5>(thrust::get<0>(tpl))),
          vt(      thrust::get<6>(thrust::get<0>(tpl)) * si::metres_per_second),
          p(       thrust::get<1>(tpl) * si::pascals),
          RH(      thrust::get<2>(tpl)),
          lambda_D(thrust::get<7>(thrust::get<0>(tpl)) * si::metres),
          lambda_K(thrust::get<8>(thrust::get<0>(tpl)) * si::metres),
          RH_max(RH_max)
        {}

        BOOST_GPU_ENABLED
        quantity<divide_typeof_helper<si::area, si::time>::type, real_t> drw2_dt(const quantity<si::area, real_t> &rw2) const
        {
          using namespace common::maxwell_mason;
          using namespace common::kappa_koehler;
          using namespace common::kelvin;
          using common::moist_air::D_0;
          using common::moist_air::K_0;
          using common::moist_air::c_pd;
          using common::transition_regime::beta;
          using common::ventil::Sh;
          using common::ventil::Nu;
#if !defined(__NVCC__)
          using std::sqrt;
#endif

          const quantity<si::length, real_t> rw  = sqrt(real_t(rw2 / si::square_metres)) * si::metres; 
          const quantity<si::volume, real_t> rw3 = rw * rw * rw;;

          // TODO: common::moist_air:: below should not be needed
          // TODO: ventilation as option
          const quantity<si::dimensionless, real_t>
            Re = common::ventil::Re(vt, rw, rhod, eta), 
            Sc = common::ventil::Sc(eta, rhod, D_0<real_t>()), // TODO? cache
            Pr = common::ventil::Pr(eta, c_pd<real_t>(), K_0<real_t>()); // TODO? cache

          const quantity<common::diffusivity, real_t> 
            D = D_0<real_t>() * beta(lambda_D / rw) * (Sh(Sc, Re) / 2);

          const quantity<common::thermal_conductivity, real_t> 
            K = K_0<real_t>() * beta(lambda_K / rw) * (Nu(Pr, Re) / 2);

          return real_t(2) * rdrdt( 
            D,
            K,
            rhod * rv, 
            T, 
            p, 
            RH > RH_max ? RH_max : RH, 
            a_w(rw3, rd3, kpa),
            klvntrm(rw, T)
          );
        }

        // backward Euler scheme:
  // rw2_new = rw2_old + f_rw2(rw2_new) * dt
  // rw2_new = rw2_old + 2 * rw * f_rw(rw2_new) * dt
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2_unitless) const
        {
          const quantity<si::area, real_t> rw2 = rw2_unitless * si::square_metres; 
          return (rw2_old + dt * drw2_dt(rw2) - rw2) / si::square_metres;
        }
      };

      template <typename real_t>
      struct advance_rw2
      {
        const real_t dt, RH_max;
        detail::config<real_t> config;

        advance_rw2(const real_t &dt, const real_t &RH_max) : dt(dt), RH_max(RH_max) {}

        BOOST_GPU_ENABLED
        real_t operator()(
          const real_t &rw2_old, 
          const thrust::tuple<thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t>, real_t, real_t> &tpl
        ) const {
#if !defined(__NVCC__)
          using std::min;
          using std::max;
          using std::cbrt;
          using std::isnan;
          using std::isinf;
#endif

          auto& tpl_in = thrust::get<0>(tpl);
          const advance_rw2_minfun<real_t> f(dt, rw2_old, tpl, RH_max); 
          const real_t drw2 = dt * f.drw2_dt(rw2_old * si::square_metres) * si::seconds / si::square_metres;

#if !defined(NDEBUG)
          if(isnan(drw2) || isinf(drw2))
          {
            // as a single command to make sure the output from parallel cuda kernels is ordered
            printf("nan/inf drw2 in cond: %g  "
              "rw2_old: %g  "
              "dt: %g  "
              "RH_max: %g  "
              "rhod: %g  "
              "rv: %g  "
              "T: %g  "
              "p: %g  "
              "RH: %g  "
              "eta: %g  "
              "rd3: %g  "
              "kpa: %g  "
              "vt: %g  "
              "lambda_D: %g  "
              "lambda_K: %g\n",
               drw2, rw2_old, dt, RH_max, 
               thrust::get<0>(tpl_in), // rhod
               thrust::get<1>(tpl_in), // rv
               thrust::get<2>(tpl_in), // T
               thrust::get<1>(tpl),    // p
               thrust::get<2>(tpl),    // RH
               thrust::get<3>(tpl_in), // eta
               thrust::get<4>(tpl_in), // rd3
               thrust::get<5>(tpl_in), // kpa
               thrust::get<6>(tpl_in), // vt
               thrust::get<7>(tpl_in), // lambda_D
               thrust::get<8>(tpl_in)  // lambda_K
            );
            assert(0);
          }
#endif

          if (drw2 == 0) return rw2_old;

          const real_t rd = cbrt(thrust::get<4>(tpl_in));
          const real_t rd2 = rd*rd;
 
          const real_t 
            a = max(rd2, rw2_old + min(real_t(0), config.cond_mlt * drw2)),
            b =          rw2_old + max(real_t(0), config.cond_mlt * drw2);

          // numerics (drw2 != 0 but a==b)
          if (a == b) return rw2_old;

          real_t fa, fb;

#if !defined(NDEBUG)
          if(a>=b)
          {
            printf("a>=b a: %g b: %g  "
              "drw2: %g  "
              "rw2_old: %g  "
              "rd2: %g  "
              "dt: %g  "
              "RH_max: %g  "
              "rhod: %g  "
              "rv: %g  "
              "T: %g  "
              "p: %g  "
              "RH: %g  "
              "eta: %g  "
              "rd3: %g  "
              "kpa: %g  "
              "vt: %g\n",
               a, b, drw2, rw2_old, rd2, dt, RH_max, thrust::get<0>(tpl_in),thrust::get<1>(tpl_in),
               thrust::get<2>(tpl_in),thrust::get<1>(tpl),thrust::get<2>(tpl),thrust::get<3>(tpl_in),
               thrust::get<4>(tpl_in),thrust::get<5>(tpl_in),thrust::get<6>(tpl_in)
            );
            assert(0);
          }
#endif

          if (drw2 > 0) 
          {
            fa = drw2; // for implicit Euler its equal to min_fun(x_old) 
            fb = f(b);
          }
          else
          {
            fa = f(a);
            fb = drw2; // for implicit Euler its equal to min_fun(x_old) 
          }

          // to store the result
          real_t rw2_new;

          // root-finding ill posed => explicit Euler 
          if (fa * fb > 0) rw2_new = rw2_old + drw2;
          // otherwise implicit Euler
          else
          {
            uintmax_t n_iter = config.n_iter;
            rw2_new = common::detail::toms748_solve(f, a, b, fa, fb, config.eps_tolerance, n_iter);
          }
          // check if it doesn't evaporate too much
          if(rw2_new < rd2) rw2_new = rd2;

#if !defined(NDEBUG)
          if(isnan(rw2_new) || isinf(rw2_new))
          {
            printf("nan/inf root in cond: %g  "
              "a: %g  "
              "b: %g  "
              "fa: %g  "
              "fb: %g\n",
              rw2_new, a, b, fa, fb
            );
            assert(0);
          }
#endif
          return rw2_new;
        }
      };
    };
  };  
};
