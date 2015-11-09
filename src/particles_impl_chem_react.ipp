// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/react.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
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
            m_SO2  = thrust::get<0>(tpl) * si::kilograms,
            m_H2O2 = thrust::get<1>(tpl) * si::kilograms,
            m_O3   = thrust::get<2>(tpl) * si::kilograms,
            m_HSO3 = thrust::get<3>(tpl) * si::kilograms,
            m_S_VI = thrust::get<4>(tpl) * si::kilograms,
            m_SO3  = thrust::get<5>(tpl) * si::kilograms,
            m_H    = thrust::get<6>(tpl) * si::kilograms;
          const quantity<si::volume, real_t> V = V_ * si::cubic_metres; 

          using namespace common::molar_mass;
          using namespace common::react;

          // TODO: optimise - do not repeat or at least do note calculate it when not needed
          // helpers for O3 reactions
          quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> 
            O3_SO2  = m_O3 / V * m_SO2  / M_SO2_H2O<real_t>()* R_S_O3_k0<real_t>(),
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
                M_SO2_H2O<real_t>() / M_O3<real_t>() * O3_SO2
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
        const thrust_device::vector<real_t> &V;
        const typename thrust_device::vector<real_t>::const_iterator &m_H; 
        const int n_part;

        // ctor
        chem_rhs(
          const thrust_device::vector<real_t> &V,
          const typename thrust_device::vector<real_t>::const_iterator &m_H
        ) :
          V(V), m_H(m_H), n_part(V.size())
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

          for (int chem_iter = chem_rhs_beg; chem_iter < chem_rhs_fin; ++chem_iter)
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
                    psi.begin() + (SO2  - chem_rhs_beg) * n_part, 
                    psi.begin() + (H2O2 - chem_rhs_beg) * n_part, 
                    psi.begin() + (O3   - chem_rhs_beg) * n_part,
                    psi.begin() + (HSO3 - chem_rhs_beg) * n_part,
                    psi.begin() + (S_VI - chem_rhs_beg) * n_part, 
                    psi.begin() + (SO3  - chem_rhs_beg) * n_part, 
                    m_H
                  )), 
                  // output
                  dot_psi.begin() + (chem_iter - chem_rhs_beg) * n_part, 
                  // op
                  chem_rhs_helper<real_t>(chem_iter)
                );
#if !defined(__NVCC__) // TODO...
                assert(std::isfinite(*thrust::min_element(
                  dot_psi.begin() + (chem_iter - chem_rhs_beg) * n_part, 
                  dot_psi.begin() + (chem_iter - chem_rhs_beg) * n_part + n_part
                )));
#endif
                break;
              default: 
                assert(false);
            }
          }
        }
      };

      template <typename real_t>
      struct chem_new_rd3
      { // recalculation of dry radii basing on created H2SO4
        const quantity<si::mass_density, real_t> chem_rho;

        // ctor
        chem_new_rd3(
          const real_t &chem_rho
        ) : 
          chem_rho(chem_rho * si::kilograms / si::cubic_metres)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(
          const thrust::tuple<real_t, real_t, real_t> &tpl
        ) 
        { 
          const quantity<si::mass, real_t>
            m_S6_old  = thrust::get<0>(tpl) * si::kilograms,     // old H2SO4
            m_S6_new  = thrust::get<1>(tpl) * si::kilograms;     // new H2SO4
          const quantity<si::volume, real_t> 
            rd3       = thrust::get<2>(tpl) * si::cubic_metres;  // old dry radii^3

          return (
            rd3 + (real_t(3./4) /
#if !defined(__NVCC__)
            pi<real_t>()
#else
            CUDART_PI
#endif
            / chem_rho) * (m_S6_new - m_S6_old)
          ) / si::cubic_metres;
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_react(const real_t &dt)
    {   
      using namespace common::molar_mass; // M-prefixed

      thrust_device::vector<real_t> &V(tmp_device_real_part);

      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      //non-equilibrium chemical reactions (oxidation)
      thrust_device::vector<real_t> &old_S_VI(tmp_device_real_part_SVI);
      // copy old H2SO4 values to allow dry radii recalculation
      thrust::copy(
        chem_bgn[S_VI], chem_end[S_VI], // from
        old_S_VI.begin()                // to
      );

      chem_stepper.do_step(
        detail::chem_rhs<real_t>(V, chem_bgn[H]), // TODO: make it an impl member field
        chem_rhs, 
        real_t(0),
        dt
      );

      // recomputing dry radii
      // TODO: using namespace for S_VI
      typedef thrust::zip_iterator<
        thrust::tuple<
          typename thrust_device::vector<real_t>::iterator, // old S_VI 
          typename thrust_device::vector<real_t>::iterator, // new_S_VI
          typename thrust_device::vector<real_t>::iterator  // rd3
        >
      > zip_it_t;

      zip_it_t 
        arg_begin(thrust::make_tuple(old_S_VI.begin(), chem_bgn[S_VI], rd3.begin())),
        arg_end(  thrust::make_tuple(old_S_VI.end(),   chem_end[S_VI], rd3.end()));
       
      thrust::transform(arg_begin, arg_end, rd3.begin(), detail::chem_new_rd3<real_t>(opts_init.chem_rho));

      //debug::print(rd3);
      //debug::print(chem_bgn[S_VI], chem_end[S_VI]);
    }
  };  
};
