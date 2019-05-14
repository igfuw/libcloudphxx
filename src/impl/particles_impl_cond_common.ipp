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
          return pow(rw2, real_t(3./2));
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

        // ctor
        BOOST_GPU_ENABLED
        advance_rw2_minfun(
          const real_t &dt,
          const real_t &rw2,
          const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t> &tpl,
          const real_t &RH_max
        ) : 
          dt(dt * si::seconds), 
          rw2_old(rw2 * si::square_metres),
          rhod(    thrust::get<0>(tpl) * si::kilograms / si::cubic_metres),
          rv(      thrust::get<1>(tpl)),
          T(       thrust::get<2>(tpl) * si::kelvins),
          p(       thrust::get<3>(tpl) * si::pascals),
          RH(      thrust::get<4>(tpl)),
          eta(     thrust::get<5>(tpl) * si::pascals * si::seconds),
          rd3(     thrust::get<6>(tpl) * si::cubic_metres),
          kpa(     thrust::get<7>(tpl)),
          vt(      thrust::get<8>(tpl) * si::metres_per_second),
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
          using common::mean_free_path::lambda_D;
          using common::mean_free_path::lambda_K;
          using common::ventil::Sh;
          using common::ventil::Nu;
          using std::sqrt;

          const quantity<si::length, real_t> rw  = sqrt(real_t(rw2 / si::square_metres)) * si::metres; 
          const quantity<si::volume, real_t> rw3 = rw * rw * rw;;

          // TODO: common::moist_air:: below should not be needed
          // TODO: ventilation as option
          const quantity<si::dimensionless, real_t>
            Re = common::ventil::Re(vt, rw, rhod, eta), 
            Sc = common::ventil::Sc(eta, rhod, D_0<real_t>()), // TODO? cache
            Pr = common::ventil::Pr(eta, c_pd<real_t>(), K_0<real_t>()); // TODO? cache

          const quantity<common::diffusivity, real_t> 
            D = D_0<real_t>() * beta(lambda_D(T)    / rw) * (Sh(Sc, Re) / 2); // TODO: cache lambdas

          const quantity<common::thermal_conductivity, real_t> 
            K = K_0<real_t>() * beta(lambda_K(T, p) / rw) * (Nu(Pr, Re) / 2);

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

        template <class tpl_tpl_t>
        BOOST_GPU_ENABLED
        real_t operator()(
          const real_t &rw2_old, 
          tpl_tpl_t tpl_tpl
        ) const {
#if !defined(__NVCC__)
          using std::min;
          using std::max;
          using std::pow;
          using std::abs;
#endif

          auto tpl = thrust::get<0>(tpl_tpl);
          //auto &tpl_revp = thrust::get<1>(tpl_tpl);

          const advance_rw2_minfun<real_t> f(dt, rw2_old, tpl, RH_max); 
          const real_t drw2 = dt * f.drw2_dt(rw2_old * si::square_metres) * si::seconds / si::square_metres;

#if !defined(NDEBUG)
          if(isnan(drw2) || isinf(drw2))
          {
            printf("nan/inf drw2 in cond: %g\n",drw2);
            printf("rw2_old: %g\n",rw2_old);
            printf("dt: %g\n",dt);
            printf("RH_max: %g\n",RH_max);
            printf("rhod: %g\n",thrust::get<0>(tpl));
            printf("rv: %g\n",thrust::get<1>(tpl));
            printf("T: %g\n",thrust::get<2>(tpl));
            printf("p: %g\n",thrust::get<3>(tpl));
            printf("RH: %g\n",thrust::get<4>(tpl));
            printf("eta: %g\n",thrust::get<5>(tpl));
            printf("rd3: %g\n",thrust::get<6>(tpl));
            printf("kpa: %g\n",thrust::get<7>(tpl));
            printf("vt: %g\n",thrust::get<8>(tpl));
            assert(0);
          }
#endif

          if (drw2 == 0) return rw2_old;

          const real_t rd2 = pow(thrust::get<6>(tpl), real_t(2./3));
 
          const real_t 
            a = max(rd2, rw2_old + min(real_t(0), config.cond_mlt * drw2)),
            b =          rw2_old + max(real_t(0), config.cond_mlt * drw2);

          // numerics (drw2 != 0 but a==b)
          if (a == b) return rw2_old;

          real_t fa, fb;

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

          // store rain evaporation rates
          if(rw2_new < rw2_old)
          {
            if(rw2_old > 4e-10)       // r_rain > 20um
              thrust::get<0>(thrust::get<1>(tpl_tpl)) += pow(rw2_old, real_t(3./2)) - pow(rw2_new, real_t(3./2));
            if(rw2_old > 6.25e-10)    // r_rain > 25um
              thrust::get<1>(thrust::get<1>(tpl_tpl)) += pow(rw2_old, real_t(3./2)) - pow(rw2_new, real_t(3./2));
            if(rw2_old > 1.024e-9)    // r_rain > 32um
              thrust::get<2>(thrust::get<1>(tpl_tpl)) += pow(rw2_old, real_t(3./2)) - pow(rw2_new, real_t(3./2));
          }

#if !defined(NDEBUG)
          if(isnan(rw2_new) || isinf(rw2_new))
          {
            printf("nan/inf root in cond: %g\n",rw2_new);
            printf("a: %g\n",a);
            printf("b: %g\n",b);
            printf("fa: %g\n",fa);
            printf("fb: %g\n",fb);
            assert(0);
          }
#endif
          return rw2_new;
        }
      };
    };
  };  
};
