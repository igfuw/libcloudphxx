// // vim:filetype=cpp
// /** @file
//   * @copyright University of Warsaw
//   * @section LICENSE
//   * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
//   */
//
// // #include <thrust/iterator/transform_iterator.h>
// #include <libcloudph++/common/maxwell-mason.hpp>
// #include <libcloudph++/common/kappa_koehler.hpp>
// #include <libcloudph++/common/kelvin_term.hpp>
// #include <libcloudph++/common/transition_regime.hpp>
// #include <libcloudph++/common/ventil.hpp>
// #include <libcloudph++/common/mean_free_path.hpp>
// #include <libcloudph++/common/detail/toms748.hpp>
//
// namespace libcloudphxx
// {
//   namespace lgrngn
//   {
//     namespace detail
//     {
//
//       template <typename real_t>
//       struct advance_ice_a_minfun
//       {
//         const quantity<si::length,            real_t> ice_a_old;
//         const quantity<si::time,              real_t> dt;
//         const quantity<si::mass_density,      real_t> rhod;
//         const quantity<si::dimensionless,     real_t> rv;
//         const quantity<si::temperature,       real_t> T;
//         const quantity<si::pressure,          real_t> p;
//         const quantity<si::dimensionless,     real_t> RH_i;
//         const quantity<si::dynamic_viscosity, real_t> eta;
//         const quantity<si::volume,            real_t> rd3;
//         const quantity<si::dimensionless,     real_t> kpa;
//         const quantity<si::velocity,          real_t> vt;
//         const quantity<si::dimensionless,     real_t> RH_max;
//         const quantity<si::length,            real_t> lambda_D;
//         const quantity<si::length,            real_t> lambda_K;
//
//         // ctor
//         BOOST_GPU_ENABLED
//         advance_ice_a_minfun(
//           const real_t &dt,
//           const real_t &ice_a,
//           const thrust::tuple<thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t>, real_t, real_t> &tpl,
//           const real_t &RH_max
//         ) :
//           dt(dt * si::seconds),
//           ice_a_old(ice_a * si::square_metres),
//           rhod(    thrust::get<0>(thrust::get<0>(tpl)) * si::kilograms / si::cubic_metres),
//           rv(      thrust::get<1>(thrust::get<0>(tpl))),
//           T(       thrust::get<2>(thrust::get<0>(tpl)) * si::kelvins),
//           eta(     thrust::get<3>(thrust::get<0>(tpl)) * si::pascals * si::seconds),
//           rd3(     thrust::get<4>(thrust::get<0>(tpl)) * si::cubic_metres),
//           kpa(     thrust::get<5>(thrust::get<0>(tpl))),
//           vt(      thrust::get<6>(thrust::get<0>(tpl)) * si::metres_per_second),
//           p(       thrust::get<1>(tpl) * si::pascals),
//           RH_i(      thrust::get<2>(tpl)),
//           lambda_D(thrust::get<7>(thrust::get<0>(tpl)) * si::metres),
//           lambda_K(thrust::get<8>(thrust::get<0>(tpl)) * si::metres),
//           RH_max(RH_max)
//         {}
//
//         BOOST_GPU_ENABLED
//         quantity<divide_typeof_helper<si::length, si::time>::type, real_t> d_ice_a_dt(const quantity<si::length, real_t> &ice_a) const
//         {
//           using namespace common::maxwell_mason;
//           using namespace common::kappa_koehler;
//           using namespace common::kelvin;
//           using common::moist_air::D_0;
//           using common::moist_air::K_0;
//           using common::moist_air::c_pd;
//           using common::transition_regime::beta;
//           using common::ventil::Sh;
//           using common::ventil::Nu;
// #if !defined(__NVCC__)
//           using std::sqrt;
// #endif
//
//           const quantity<si::dimensionless, real_t>
//             Re = common::ventil::Re(vt, ice_a, rhod, eta),
//             Sc = common::ventil::Sc(eta, rhod, D_0<real_t>()), // TODO? cache
//             Pr = common::ventil::Pr(eta, c_pd<real_t>(), K_0<real_t>()); // TODO? cache
//
//           const quantity<common::diffusivity, real_t>
//             D = D_0<real_t>() * beta(lambda_D / ice_a) * (Sh(Sc, Re) / 2);
//
//           const quantity<common::thermal_conductivity, real_t>
//             K = K_0<real_t>() * beta(lambda_K / ice_a) * (Nu(Pr, Re) / 2);
//
//           return da_dt(
//             D,
//             K,
//             rhod * rv,
//             T,
//             p,
//             RH_i > RH_max ? RH_max : RH_i
//           );
//         }
//
//         BOOST_GPU_ENABLED
//         real_t operator()(const real_t &ice_a_unitless) const
//         {
//           const quantity<si::length, real_t> ice_a = ice_a_unitless * si::metres;
//           return (ice_a_old + dt * d_ice_a_dt(ice_a) - ice_a) / si::metres;
//         }
//       };
//
//
//
//       template <typename real_t>
//       struct advance_ice_axis
//       {
//         const real_t dt, RH_max;
//
//         advance_ice_axis(const real_t &dt, const real_t &RH_max)
//           : dt(dt), RH_max(RH_max) {}
//
//         BOOST_GPU_ENABLED
//         thrust::tuple<real_t, real_t> operator()(
//             const thrust::tuple<real_t, real_t> &ac_old,
//             const thrust::tuple<
//                 thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t>,
//                 real_t, real_t> &tpl
//         ) const
//         {
// #if !defined(__NVCC__)
//           using std::max;
//           using std::isnan;
//           using std::isinf;
// #endif
//           const real_t a_old = thrust::get<0>(ac_old);
//           const real_t c_old = thrust::get<1>(ac_old);
//
//           // Skip liquid droplets
//           if (a_old <= real_t(0) || c_old <= real_t(0))
//             return ac_old;
//
//           advance_ice_a_minfun<real_t> f_a(dt, a_old * a_old, tpl, RH_max);
//           advance_ice_c_minfun<real_t> f_c(dt, c_old * c_old, tpl, RH_max);
//
//           const real_t da_dt = (f_a.drw2_dt(a_old * a_old * si::square_metres) / (2 * a_old * si::metres))
//                                 * si::seconds / si::metres;
//           const real_t dc_dt = (f_c.drw2_dt(c_old * c_old * si::square_metres) / (2 * c_old * si::metres))
//                                 * si::seconds / si::metres;
//
//           // to store the result
//           real_t a_new, c_new;
//
//           if (da_dt == 0 && dc_dt == 0) return ac_old;
//
//           const real_t rd = cbrt(thrust::get<4>(tpl_in));
//
//           const real_t
//             a = max(rd2, rw2_old + min(real_t(0), cond_mlt * drw2)),
//             b =          rw2_old + max(real_t(0), cond_mlt * drw2);
//
//           // numerics (drw2 != 0 but a==b)
//           if (a == b)
//           {
//             if constexpr (apply)
//               return rw2_old;
//             else
//               return real_t(0);
//           }
//
//           real_t fa, fb;
//
//           if (drw2 > 0)
//           {
//             fa = drw2; // for implicit Euler its equal to min_fun(x_old)
//             fb = f(b);
//           }
//           else
//           {
//             fa = f(a);
//             fb = drw2; // for implicit Euler its equal to min_fun(x_old)
//           }
//
//           // root-finding ill posed => explicit Euler
//           if (fa * fb > 0) rw2_new = rw2_old + drw2;
//           // otherwise implicit Euler
//           else
//           {
//             auto _n_iter = n_iter; // we need a copy because toms748_solve expects non-const ref (n_iter is modified by it)
//             rw2_new = common::detail::toms748_solve(f, a, b, fa, fb, eps_tolerance, _n_iter);
//           }
//
//           // check if it doesn't evaporate too much
//           if(rw2_new < rd2) rw2_new = rd2;
//
//           return thrust::make_tuple(a_new, c_new, rho_new);
//         }
//       };
//
//     };
//   };
// };
