// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/detail/toms748.hpp>
#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/henry.hpp>
#include <libcloudph++/common/dissoc.hpp>
#include <libcloudph++/common/react.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct chem_vol_fun
      { // calculate drop volume
        const real_t pi;

        // ctor (pi() is not a __device__ function...)
        chem_vol_fun() :
          pi(boost::math::constants::pi<real_t>())
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2) const
        {
#if !defined(__NVCC__)
	  using std::pow;
#endif
          return real_t(4./3) * pi * (pow(rw2, real_t(3./2)));
        }
      };

      template <typename real_t>
      struct chem_Henry_fun
      { // gas absorption into cloud droplets (Henrys law)
        const quantity<common::amount_over_volume_over_pressure, real_t> H;
        const quantity<common::mass_over_amount, real_t> M;
        const quantity<si::dimensionless, real_t> c;

        // ctor
        chem_Henry_fun(
          const quantity<common::amount_over_volume_over_pressure, real_t> &H,
          const quantity<common::mass_over_amount, real_t> &M,
          const quantity<si::dimensionless, real_t> &c
        ) : 
          H(H), M(M), c(c)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &V, const real_t &p) const
        { 
          return (
            H                       // Henry's law
            * c * (p * si::pascals) // volume concentration -> partial pressure
            * V * si::cubic_metres  // drop volume 
	    * M                     // moles -> kilograms
          ) / si::kilograms;
        }
      };

      template <typename real_t>
      struct chem_minfun
      { // helper function for maintaining electroneutrality
        // returns the difference between the mass of H+ ions already present 
        // and the mass of H+ ions needed to be added (due to dissociation)
        // in order to maintain electroneutrality
	const quantity<si::mass, real_t> &m_SO2, &m_S_VI;
	const quantity<si::volume, real_t> &V;
        
        // ctor
        BOOST_GPU_ENABLED
        chem_minfun(
	  const quantity<si::mass, real_t> &m_SO2,
	  const quantity<si::mass, real_t> &m_S_VI,
	  const quantity<si::volume, real_t> &V
        ) :
          m_SO2(m_SO2), m_S_VI(m_S_VI), V(V)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &arg) const
        {
	  using namespace common::molar_mass;
	  using namespace common::dissoc;

          const quantity<si::mass, real_t> m_H = arg * si::kilograms;

          return (-m_H + M_H<real_t>() * (
            // HSO3 to SO3 dissoctation 
            real_t(2) * // "2-" ion //TODO is it correct if we already start from HSO3- ion?
            K_SO2<real_t>() * K_HSO3<real_t>() / M_SO2<real_t>() * m_SO2 
            * (M_H<real_t>() * M_H<real_t>()) 
            * (V * V) 
            / (m_H * m_H)
            +
            // H2O*SO2 to HSO3 dissociation
            K_SO2<real_t>() / M_SO2<real_t>() * m_SO2 * M_H<real_t>() * V / m_H
            +
            // dissociation of pure water 
            K_H2O<real_t>() * M_H<real_t>() * (V*V) / m_H
            +
            // dissociation of S_VI to HSO4
            m_H * m_S_VI / V / M_H2SO4<real_t>() / M_H<real_t>()
	    / ((m_H) / M_H<real_t>() / V + K_HSO4<real_t>())
            +
            // dissociation of HSO4 to SO4
            real_t(2) * // "2-" ion
            K_HSO4<real_t>() * m_S_VI 
            / M_H2SO4<real_t>()
            / ((m_H) / M_H<real_t>() / V + K_HSO4<real_t>())


          )) / si::kilograms;
        }
      };

      template <typename real_t>
      struct chem_electroneutral // TODO: does it have to be a struct/functor - perhaps ordinary function would suffice?
      { // uses toms748 scheme to solve for mass of H+ after dissociation
        // that satisfies electroneutrality
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) const
        {
          using namespace common::molar_mass;

          const quantity<si::mass, real_t>
            m_SO2  = thrust::get<0>(tpl) * si::kilograms,
            m_S_VI = thrust::get<1>(tpl) * si::kilograms; 
          const quantity<si::volume, real_t> 
            V      = thrust::get<2>(tpl) * si::cubic_metres;
          
          // left side for search in toms748
          real_t m_H_pure = ((real_t(1e-7 * 1e3) * si::moles / si::cubic_metres) * V * M_H<real_t>()) / si::kilograms;

          real_t m_H = common::detail::toms748_solve(
	    detail::chem_minfun<real_t>(m_SO2, m_S_VI, V),
            m_H_pure, // min -> (pure water)
	    real_t(1e-10) // max -> TODO
	  ); 
          //std::cerr << "  " << m_H_pure << " ... " << m_H << " ... " << "TODO" << std::endl;
          // TODO: asserts for K = f(m_H, m_...)
          return m_H;
        }
      };
 
      template <typename real_t, int chem_iter>
      struct chem_dissoc_diag // TODO: does it have to be a struct/functor - perhaps ordinary function would suffice?
      { // basing on the mass of H+ ions from chem_electroneutral
        // it recalculatess the mass of all other chemical species after dissociation
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::volume, real_t> 
            V      = thrust::get<0>(tpl) * si::cubic_metres;
          const quantity<si::mass, real_t>
            m_H    = thrust::get<1>(tpl) * si::kilograms,
            m_SO2  = thrust::get<2>(tpl) * si::kilograms,
            m_HSO3 = thrust::get<3>(tpl) * si::kilograms,
            m_S_VI = thrust::get<4>(tpl) * si::kilograms;

          switch (chem_iter)
          {
	    using namespace common::dissoc;     // K-prefixed
	    using namespace common::molar_mass; // M-prefixed

            case OH:
              return (
                M_H<real_t>() * M_OH<real_t>() 
                * K_H2O<real_t>()               // note: dissociation constant for pure water 
                * V * V / m_H                   // is actually k*[H2O] (Seinfeld and Pandis p 345)
              ) / si::kilograms; 
            case HSO3:
              return (
                V / m_H * m_SO2
                * M_HSO3<real_t>() * K_SO2<real_t>() * M_H<real_t>() / M_SO2<real_t>()
              ) / si::kilograms;
            case SO3:
              return (
                V
                * M_SO3<real_t>()
                * K_HSO3<real_t>() 
                * m_HSO3 / m_H
                * M_H<real_t>() / M_HSO3<real_t>()
              ) / si::kilograms;
            case HSO4:
              return (
                M_HSO4<real_t>() / M_H<real_t>() / M_H2SO4<real_t>()
                * m_H
                * m_S_VI
                / V 
                / ( m_H / M_H<real_t>() / V + K_HSO4<real_t>())
              ) / si::kilograms;
            case SO4:
              return (
                M_SO4<real_t>() / M_H2SO4<real_t>()
                * K_HSO4<real_t>() 
                * m_S_VI
                / (m_H / M_H<real_t>() / V + K_HSO4<real_t>())
              ) / si::kilograms;
            default:
              assert(false);
              return 0;
          }
        }
      };

      template <typename real_t>
      struct chem_rhs_helper
      { // returns the change in mass per second of each chemical compounds
        // due to oxidation reaction
        const int chem_iter;
        
        // ctor
        chem_rhs_helper(const int &chem_iter) : 
          chem_iter(chem_iter)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(
          const real_t &V_, 
          const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t> &tpl
        ) const
        {
          const quantity<si::mass, real_t>
            m_SO2  = thrust::get<SO2 >(tpl) * si::kilograms,
            m_H2O2 = thrust::get<H2O2>(tpl) * si::kilograms,
            m_O3   = thrust::get<O3  >(tpl) * si::kilograms,
            m_HSO3 = thrust::get<HSO3>(tpl) * si::kilograms,
            m_S_VI = thrust::get<S_VI>(tpl) * si::kilograms,
            m_SO3  = thrust::get<SO3 >(tpl) * si::kilograms,
            m_H    = thrust::get<H   >(tpl) * si::kilograms;
          const quantity<si::volume, real_t> V = V_ * si::cubic_metres; 

          using namespace common::molar_mass;
          using namespace common::react;

          // TODO: optimise - do not repeat or at least do note calculate it when not needed
          // helpers for O3 reactions
          quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> 
            O3_SO2  = m_O3 / V * m_SO2  / M_SO2<real_t>()  * R_S_O3_k0<real_t>(),
            O3_HSO3 = m_O3 / V * m_HSO3 / M_HSO3<real_t>() * R_S_O3_k1<real_t>(),
            O3_SO3  = m_O3 / V * m_SO3  / M_SO3<real_t>()  * R_S_O3_k2<real_t>();


	  // helper for H2O2 reactions
	  quantity<divide_typeof_helper<si::amount, si::time>::type, real_t> 
            H2O2_HSO3 = R_S_H2O2_k<real_t>() / (V*V) 
	      * m_H2O2 / M_H2O2<real_t>()
	      * m_H    / M_H<real_t>()
	      * m_HSO3 / M_HSO3<real_t>()
	      / (real_t(1) + R_S_H2O2_K<real_t>() * m_H / M_H<real_t>() / V);

          switch (chem_iter)
          {
            case S_VI:
              return (
                M_H2SO4<real_t>() / M_O3<real_t>() * (O3_SO2 + O3_HSO3 + O3_SO3)
                +
                M_H2SO4<real_t>() * H2O2_HSO3
              ) / si::kilograms * si::seconds;
            case H2O2:
              return -(
                M_H2O2<real_t>() * H2O2_HSO3
              ) / si::kilograms * si::seconds;
            case O3:
              return -(
                O3_SO2 + O3_HSO3 + O3_SO3
              ) / si::kilograms * si::seconds;
            case SO2:
              return -(
                M_SO2<real_t>() / M_O3<real_t>() * O3_SO2
              ) / si::kilograms * si::seconds;
            case HSO3:
              return -(
		M_HSO3<real_t>() / M_O3<real_t>() * O3_HSO3
		+
		M_HSO3<real_t>() * H2O2_HSO3
              ) / si::kilograms * si::seconds;
            case SO3:
              return (
                - M_SO3<real_t>() / M_O3<real_t>() * O3_SO3
              ) / si::kilograms * si::seconds;
            default:
              assert(false);
              return 0;
          }
        }
      };

      // functor called by odeint
      template <typename real_t>
      struct chem_rhs
      {
        const thrust_device::vector<real_t> &V, &chem_equil;
        const int n_part;

        // ctor
        chem_rhs(
          const thrust_device::vector<real_t> &V,
          const thrust_device::vector<real_t> &chem_equil
        ) :
          V(V), chem_equil(chem_equil), n_part(V.size())
        {}

        void operator()(
          const thrust_device::vector<real_t> &psi, 
          thrust_device::vector<real_t> &dot_psi,
          const real_t /* t */
        )
        {
          assert(dot_psi.size() == psi.size());
          typedef thrust::zip_iterator<
            thrust::tuple<
              // those in psi...
              typename thrust_device::vector<real_t>::const_iterator, // SO2
              typename thrust_device::vector<real_t>::const_iterator, // H2O2
              typename thrust_device::vector<real_t>::const_iterator, // O3
              typename thrust_device::vector<real_t>::const_iterator, // HSO3
              typename thrust_device::vector<real_t>::const_iterator, // S_VI
              typename thrust_device::vector<real_t>::const_iterator, // SO3
              // ... but not only!
              typename thrust_device::vector<real_t>::const_iterator  // H
            >
          > zip_it_t;

          for (int chem_iter = 0; chem_iter < chem_rhs_n; ++chem_iter)
          {
            switch (chem_iter)
            {
              case S_VI: 
              case H2O2:
              case O3:
              case SO2:
              case HSO3:
              case SO3:
                thrust::transform(
                  // input - 1st arg
                  V.begin(), V.end(),                
                  // input - 2nd arg
	          zip_it_t(thrust::make_tuple(
                    psi.begin() + 0 * n_part, 
                    psi.begin() + 1 * n_part, 
                    psi.begin() + 2 * n_part,
                    psi.begin() + 3 * n_part,
                    psi.begin() + 4 * n_part, 
                    psi.begin() + 5 * n_part, 
                    chem_equil.begin() + 0 * n_part // H - hardcoded! has to comply wiith enum definition :(
                  )), 
                  // output
                  dot_psi.begin() + chem_iter * n_part, 
                  // op
                  chem_rhs_helper<real_t>(chem_iter)
                );
                break;
              default: 
                assert(false);
            }
          }
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem(
      const real_t &dt,
      const std::vector<real_t> &chem_gas,
      const bool &chem_dsl, const bool &chem_dsc, const bool &chem_rct
    )
    {   
      using namespace common::henry;      // H-prefixed
      using namespace common::molar_mass; // M-prefixed

      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      // 0/4: calculating drop volumes
      thrust_device::vector<real_t> &V(tmp_device_real_part);
      thrust::transform(
        rw2.begin(), rw2.end(),         // input
        V.begin(),                      // output 
        detail::chem_vol_fun<real_t>()   // op
      );

      if (chem_dsl == true){  //TODO move to a separate function and then move the logic to opts particle_step 
        // 1/4: equilibrium stuff: gas absortption
        // TODO: open/close system logic
        // TODO: K=K(T)
        assert(SO2 == 0 && H2O2 == 1 && O3 == 2);
        static const quantity<common::amount_over_volume_over_pressure, real_t> H_[chem_gas_n] = {
	  H_SO2<real_t>(), H_H2O2<real_t>(), H_O3<real_t>()
        };
        static const quantity<common::mass_over_amount, real_t> M_[chem_gas_n] = {
          M_SO2<real_t>(), M_H2O2<real_t>(), M_O3<real_t>()
        };
        for (int i = 0; i < chem_gas_n; ++i)
        {
	  thrust::transform(
	    V.begin(), V.end(),                                        // input - 1st arg
	    thrust::make_permutation_iterator(p.begin(), ijk.begin()), // input - 2nd arg
	    chem_bgn[i],                                               // output
	    detail::chem_Henry_fun<real_t>(H_[i], M_[i], chem_gas[i])     // op
	  );
        }
      }

      if (chem_dsc == true){  //TODO move to a separate function and then move the logic to opts particle_step 
        // 2/4: equilibrium stuff: dissociation
        { // H+ ions after dissociation so that drops remain electroneutral
          typedef thrust::zip_iterator<
            thrust::tuple<
              typename thrust_device::vector<real_t>::iterator, // SO2
              typename thrust_device::vector<real_t>::iterator, // S_VI
              typename thrust_device::vector<real_t>::iterator  // V
            >
          > zip_it_t;

          thrust::transform(
	    zip_it_t(thrust::make_tuple(chem_bgn[SO2], chem_bgn[S_VI], V.begin())),  // input - begin
	    zip_it_t(thrust::make_tuple(chem_end[SO2], chem_end[S_VI], V.end())  ),  // input - end
	    chem_bgn[H],                                                             // output
	    detail::chem_electroneutral<real_t>()                                          // op
	  );
        }
        { // diagnose the rest of ions basing on H+
          typedef thrust::zip_iterator<
            thrust::tuple<
              typename thrust_device::vector<real_t>::iterator, // V
              typename thrust_device::vector<real_t>::iterator, // H 
              typename thrust_device::vector<real_t>::iterator, // SO2
              typename thrust_device::vector<real_t>::iterator, // HSO3
              typename thrust_device::vector<real_t>::iterator  // S_VI
            >
          > zip_it_t;

          zip_it_t 
            arg_begin(thrust::make_tuple(V.begin(), chem_bgn[H], chem_bgn[SO2], chem_bgn[HSO3], chem_bgn[S_VI])),
            arg_end(  thrust::make_tuple(V.end(),   chem_end[H], chem_end[SO2], chem_end[HSO3], chem_end[S_VI]));
        
          thrust::transform(arg_begin, arg_end, chem_bgn[OH  ], detail::chem_dissoc_diag<real_t, OH  >());
          thrust::transform(arg_begin, arg_end, chem_bgn[HSO3], detail::chem_dissoc_diag<real_t, HSO3>()); // note: has to be computed before SO3
          thrust::transform(arg_begin, arg_end, chem_bgn[SO3 ], detail::chem_dissoc_diag<real_t, SO3 >());
          thrust::transform(arg_begin, arg_end, chem_bgn[HSO4], detail::chem_dissoc_diag<real_t, HSO4>());
          thrust::transform(arg_begin, arg_end, chem_bgn[SO4 ], detail::chem_dissoc_diag<real_t, SO4 >());
        }
      }

      if (chem_rct == true){  //TODO move to a separate function and then move the logic to opts particle_step 
        // 3/4: non-equilibrium stuff
        {
          chem_stepper.do_step(
            detail::chem_rhs<real_t>(V, chem_equil), // TODO: make it an impl member field
            chem_noneq, 
            real_t(0),
            dt
          );
        }
      }

      // 4/4: recomputing dry radii
      {
        namespace arg = thrust::placeholders;
        // TODO: using namespace for S_VI
        thrust::transform(
          chem_bgn[S_VI], chem_end[S_VI],                               // input
          rd3.begin(),                                                  // output
          (real_t(3./4) / pi<real_t>() / opts_init.chem_rho) * arg::_1  // op
        );

        //debug::print(rd3)
        //debug::print(chem_bgn[S_VI], chem_end[S_VI])
      };
    }
  };  
};
