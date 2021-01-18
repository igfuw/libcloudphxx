// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/react.hpp>
#include <libcloudph++/common/dissoc.hpp>

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
        const real_t dt;
        
        // ctor
        chem_rhs_helper(const int &chem_iter, const real_t dt) : 
          chem_iter(chem_iter), dt(dt)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(
          const real_t &V_, 
          const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t> &tpl
        ) const
        {
          const quantity<si::volume, real_t>      V = V_ * si::cubic_metres; 
          const quantity<si::temperature, real_t> T = thrust::get<0>(tpl) * si::kelvins;

          const quantity<si::mass, real_t>
            m_S_IV = thrust::get<1>(tpl) * si::kilograms,
            m_S_VI = thrust::get<2>(tpl) * si::kilograms,
            m_H2O2 = thrust::get<3>(tpl) * si::kilograms,
            m_O3   = thrust::get<4>(tpl) * si::kilograms,
            m_H    = thrust::get<5>(tpl) * si::kilograms;

          using namespace common::molar_mass;
          using namespace common::dissoc;
          using namespace common::react;

          // helper for H+ concentration
          quantity<common::amount_over_volume, real_t> conc_H;
          conc_H = m_H / M_H<real_t>() / V;

          //helpers for dissociation (temperature dependance)
          quantity<common::amount_over_volume, real_t> Kt_SO2, Kt_HSO3;
          Kt_SO2  = K_temp(T, K_SO2<real_t>(),  dKR_SO2<real_t>());
          Kt_HSO3 = K_temp(T, K_HSO3<real_t>(), dKR_HSO3<real_t>());

          // helpers for reactions (temperature dependance)
          quantity<common::volume_over_amount_over_time, real_t> R_O3_k0, R_O3_k1, R_O3_k2;
          quantity<common::volume_square_over_amount_square_over_time, real_t> R_H2O2_k;
          R_O3_k0  = R_temp_O3(T, R_S_O3_k0<real_t>(),  dER_O3_k0<real_t>());
          R_O3_k1  = R_temp_O3(T, R_S_O3_k1<real_t>(),  dER_O3_k1<real_t>());
          R_O3_k2  = R_temp_O3(T, R_S_O3_k2<real_t>(),  dER_O3_k2<real_t>());
          R_H2O2_k = R_temp_H2O2(T, R_S_H2O2_k<real_t>(), dER_H2O2_k<real_t>()); 

          // helper for O3 reaction
          quantity<divide_typeof_helper<si::amount, si::time>::type, real_t> 
            O3_react = V * m_O3 / M_O3<real_t>() / V * m_S_IV / M_SO2_H2O<real_t>() / V
                       / (real_t(1) + Kt_SO2 / conc_H + Kt_SO2 * Kt_HSO3 / conc_H / conc_H)
                       * (R_O3_k0 + R_O3_k1 * Kt_SO2 / conc_H + R_O3_k2 * Kt_SO2 * Kt_HSO3 / conc_H / conc_H);
            
            // check if reaction rate won't take more O3 than there is
            O3_react = (O3_react * dt * si::seconds < m_O3 / M_O3<real_t>()) ? O3_react : 
                         m_O3 / M_O3<real_t>() / (dt * si::seconds);
            // check if reaction rate won't take more S_IV than there is
            O3_react = (O3_react * dt * si::seconds < m_S_IV / M_SO2_H2O<real_t>()) ? O3_react : 
                         m_S_IV / M_SO2_H2O<real_t>() / (dt * si::seconds);
 
	  // helper for H2O2 reaction
	  quantity<divide_typeof_helper<si::amount, si::time>::type, real_t>
            H2O2_react = V * R_H2O2_k * Kt_SO2 
                         * m_H2O2 / M_H2O2<real_t>() / V
                         * m_S_IV / M_SO2_H2O<real_t>() / V
                         / (real_t(1) + Kt_SO2 / conc_H + Kt_SO2 * Kt_HSO3 / conc_H / conc_H)
                         / (real_t(1) + R_S_H2O2_K<real_t>() * conc_H);

            // check if reaction rate won't take more H2O2 than there is
            H2O2_react = (H2O2_react * dt * si::seconds < m_H2O2 / M_H2O2<real_t>()) ? H2O2_react : 
                           m_H2O2 / M_H2O2<real_t>() / (dt * si::seconds);
            // check if reaction rate won't take more S_IV than there is (this silently gives precedence to O3 reaction)
            H2O2_react = (H2O2_react * dt * si::seconds < m_S_IV / M_SO2_H2O<real_t>() - O3_react * dt * si::seconds) ? 
                          H2O2_react : m_S_IV / M_SO2_H2O<real_t>() / (dt * si::seconds) - O3_react;

          switch (chem_iter)
          {
            case SO2:
              return -(
                M_SO2_H2O<real_t>() * (O3_react + H2O2_react)
              ) / si::kilograms * si::seconds;
            case S_VI:
              return (
                M_H2SO4<real_t>() * (O3_react + H2O2_react)
              ) / si::kilograms * si::seconds;
            case H2O2:
              return -(
                M_H2O2<real_t>() * H2O2_react
              ) / si::kilograms * si::seconds;
            case O3:
              return -(
                M_O3<real_t>() * O3_react
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
        const real_t dt;
        const thrust_device::vector<real_t> &V;
        const thrust_device::vector<unsigned int> &chem_flag;
 
        typedef thrust::permutation_iterator<
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<thrust_size_t>::iterator
        > pi_t;
        const pi_t &T;
       
        const typename thrust_device::vector<real_t>::const_iterator &m_H; 
        const int n_part;

        // ctor
        chem_rhs(
          const real_t dt,
          const thrust_device::vector<real_t> &V,
          const pi_t &T,
          const typename thrust_device::vector<real_t>::const_iterator &m_H,
          const thrust_device::vector<unsigned int> &chem_flag
        ) :
          dt(dt), V(V), T(T), m_H(m_H), n_part(V.size()), chem_flag(chem_flag)
        {}

        void operator()(
          const thrust_device::vector<real_t> &psi, 
          thrust_device::vector<real_t> &dot_psi,
          const real_t /* t */
        )
        {
          thrust::fill(dot_psi.begin(), dot_psi.end(), real_t(0));
          assert(dot_psi.size() == psi.size());

          typedef thrust::zip_iterator<
            thrust::tuple<
              pi_t, // T
              // those in psi...
              typename thrust_device::vector<real_t>::const_iterator, // S_IV
              typename thrust_device::vector<real_t>::const_iterator, // S_VI
              typename thrust_device::vector<real_t>::const_iterator, // H2O2
              typename thrust_device::vector<real_t>::const_iterator, // O3
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
                thrust::transform_if(
                  // input - 1st arg
                  V.begin(), V.end(),                
                  // input - 2nd arg
	          zip_it_t(thrust::make_tuple(
                    T,
                    psi.begin() + (SO2  - chem_rhs_beg) * n_part, 
                    psi.begin() + (S_VI - chem_rhs_beg) * n_part, 
                    psi.begin() + (H2O2 - chem_rhs_beg) * n_part, 
                    psi.begin() + (O3   - chem_rhs_beg) * n_part,
                    m_H
                  )), 
                  // chemical reactions are only done for selected droplets 
                  // (with wet radius significantly bigger than dry radius)
                  // to avoid problems in activation when dry radius (due to chemistry) 
                  // is bigger than predicted wet radius in condensation
                  // 
                  // stencil 
                  chem_flag.begin(),
                  // output
                  dot_psi.begin() + (chem_iter - chem_rhs_beg) * n_part, 
                  // op
                  chem_rhs_helper<real_t>(chem_iter, dt),
                  // condition
                  thrust::identity<unsigned int>()
                );

#if !defined(__NVCC__)
                using boost::math::isfinite;
#endif
                assert(isfinite(*thrust::min_element(
                  dot_psi.begin() + (chem_iter - chem_rhs_beg) * n_part, 
                  dot_psi.begin() + (chem_iter - chem_rhs_beg) * n_part + n_part
                )));
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
        const real_t chem_rho;

        // ctor
        chem_new_rd3(
          const real_t &chem_rho
        ) : 
          chem_rho(chem_rho)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(
          const thrust::tuple<real_t, real_t, real_t> &tpl
        ) 
        { 
          //const quantity<si::mass, real_t>
          const real_t 
            m_S6_old  = thrust::get<0>(tpl),// * si::kilograms,     // old H2SO4
            m_S6_new  = thrust::get<1>(tpl);// * si::kilograms;     // new H2SO4
          //const quantity<si::volume, real_t>
          const real_t  
            rd3       = thrust::get<2>(tpl);// * si::cubic_metres;  // old dry radii^3

          return 
            rd3 + (real_t(3./4) /
#if !defined(__NVCC__)
            pi<real_t>()
#else
            CUDART_PI
#endif
            / chem_rho) * (m_S6_new - m_S6_old);
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_react(const real_t &dt)
    {   
      using namespace common::molar_mass; // M-prefixed

      thrust_device::vector<real_t> &V(tmp_device_real_part);
      thrust_device::vector<unsigned int> &chem_flag(tmp_device_n_part);

      //non-equilibrium chemical reactions (oxidation)
      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      thrust_device::vector<real_t> &old_S_VI(tmp_device_real_part1);

      // copy old H2SO4 values to allow dry radii recalculation
      thrust::copy(
        chem_bgn[S_VI], chem_end[S_VI], // from
        old_S_VI.begin()                // to
      );

      // do chemical reactions
      chem_stepper.do_step(
        detail::chem_rhs<real_t>(
          dt,
          V,
          thrust::make_permutation_iterator(T.begin(), ijk.begin()), 
          chem_bgn[H], 
          chem_flag
        ), // TODO: make it an impl member field
        chem_rhs, 
        real_t(0),
        dt
      );

      assert(opts_init.chem_rho != 0);

      // recompute dry radii
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
 
      // do oxidation reaction only for droplets that are "big enough" (marked by chem_flag) 
      thrust::transform_if(
        arg_begin, arg_end,                                //input first arg 
        chem_flag.begin(),                                 //stencil
        rd3.begin(),                                       //output
        detail::chem_new_rd3<real_t>(opts_init.chem_rho),  //op
        thrust::identity<unsigned int>()                   //condition
      );

#if !defined(__NVCC__)
      using boost::math::isfinite;
#endif
      for (int i = 0; i < chem_gas_n; ++i){
        assert(isfinite(*thrust::min_element(chem_bgn[i], chem_end[i])));
      }

      assert(isfinite(*thrust::min_element(rd3.begin(), rd3.end())));
    }
  };  
};
