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
        const quantity<si::temperature, real_t> dHR;
        const quantity<common::mass_over_amount, real_t> M;
        const quantity<si::dimensionless, real_t> c;
        const quantity<common::diffusivity, real_t> D;
        const quantity<si::dimensionless, real_t> acc_coeff;
        const quantity<si::time, real_t> dt;

        // ctor
        chem_Henry_fun(
          const quantity<common::amount_over_volume_over_pressure, real_t> &H,
          const quantity<si::temperature, real_t> &dHR,
          const quantity<common::mass_over_amount, real_t> &M,
          const quantity<si::dimensionless, real_t> &c,
          const quantity<common::diffusivity, real_t> &D,
          const quantity<si::dimensionless, real_t> &acc_coeff,
          const quantity<si::time, real_t> &dt
        ) : 
          H(H), dHR(dHR), M(M), c(c), D(D), acc_coeff(acc_coeff), dt(dt)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &V, const thrust::tuple<real_t, real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::pressure, real_t>    p      = thrust::get<0>(tpl) * si::pascals; 
          const quantity<si::temperature, real_t> T      = thrust::get<1>(tpl) * si::kelvins;     
          const quantity<si::mass, real_t>        m_old  = thrust::get<2>(tpl) * si::kilograms;     
          const quantity<si::area, real_t>        rw2    = thrust::get<3>(tpl) * si::metres * si::metres;  

          //implicit solution to the eq. 8.22 from chapter 8.4.2 in Peter Warneck Chemistry of the Natural Atmosphere  
          return (
            (m_old + M * V * si::cubic_metres * c * p / T / common::moist_air::kaBoNA<real_t>() 
                     * dt * common::henry::mass_trans(rw2, D, acc_coeff, T, M)) 
            /
            (real_t(1.) 
              + common::henry::mass_trans(rw2, D, acc_coeff, T, M) * dt 
                / common::henry::H_temp(T, H, dHR) / common::moist_air::kaBoNA<real_t>() / T
            )
          ) / si::kilograms;
        }
      };

      template <typename real_t>
      struct chem_save_moles_2
      {
        const quantity<common::mass_over_amount, real_t> M1, M2, M3; //molar mass
	
        // ctor
        BOOST_GPU_ENABLED
        chem_save_moles_2(
          const quantity<common::mass_over_amount, real_t> &M1,
          const quantity<common::mass_over_amount, real_t> &M2,
          const quantity<common::mass_over_amount, real_t> &M3
        ) :
          M1(M1), M2(M2), M3(M3)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::mass, real_t> m1  = thrust::get<0>(tpl) * si::kilograms;     
          const quantity<si::mass, real_t> m2  = thrust::get<1>(tpl) * si::kilograms;     
          const quantity<si::mass, real_t> m3  = thrust::get<2>(tpl) * si::kilograms;
          
          return (m1 / M1 + m2 / M2 + m3 / M3) / si::moles;
        }
      };

      template <typename real_t>
      struct chem_save_moles_1
      {
        const quantity<common::mass_over_amount, real_t> M1, M2; //molar mass
	
        // ctor
        BOOST_GPU_ENABLED
        chem_save_moles_1(
          const quantity<common::mass_over_amount, real_t> &M1,
          const quantity<common::mass_over_amount, real_t> &M2
        ) :
          M1(M1), M2(M2)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t> &tpl) const
        {
          const quantity<si::mass, real_t> m1  = thrust::get<0>(tpl) * si::kilograms;     
          const quantity<si::mass, real_t> m2  = thrust::get<1>(tpl) * si::kilograms;     
          
          return (m1 / M1 + m2 / M2) / si::moles;
        }
      };
 
      template <typename real_t>
      struct chem_minfun
      { // helper function for maintaining electroneutrality
        // returns the difference between the mass of H+ ions already present 
        // and the mass of H+ ions needed to be added (due to dissociation)
        // in order to maintain electroneutrality
	const quantity<si::amount, real_t> n_SO2_old, n_CO2_old, n_HNO3_old, n_NH3_old;
	const quantity<si::mass, real_t> m_S_VI;
	const quantity<si::volume, real_t> V;
        
        // ctor
        BOOST_GPU_ENABLED
        chem_minfun(
	  const quantity<si::amount, real_t> &n_SO2_old,
	  const quantity<si::amount, real_t> &n_CO2_old,
	  const quantity<si::amount, real_t> &n_HNO3_old,
	  const quantity<si::amount, real_t> &n_NH3_old,
	  const quantity<si::mass, real_t> &m_S_VI,
	  const quantity<si::volume, real_t> &V
        ) :
          n_SO2_old(n_SO2_old), n_CO2_old(n_CO2_old), n_HNO3_old(n_HNO3_old), n_NH3_old(n_NH3_old), m_S_VI(m_S_VI), V(V)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &arg) const
        {
	  using namespace common::molar_mass;
	  using namespace common::dissoc;

          const quantity<si::mass, real_t> m_H = arg * si::kilograms;

          quantity<common::amount_over_volume, real_t> conc_H, CO2_helper, SO2_helper;

          conc_H = m_H / M_H<real_t>() / V;

          CO2_helper = n_CO2_old / V 
                       / (real_t(1) + K_CO2<real_t>() / conc_H + K_CO2<real_t>() * K_HCO3<real_t>() / conc_H / conc_H);
          SO2_helper = n_SO2_old / V 
                       / (real_t(1) + K_SO2<real_t>() / conc_H + K_SO2<real_t>() * K_HSO3<real_t>() / conc_H / conc_H);

          return (-m_H + M_H<real_t>() * (
            // dissociation of pure water 
            K_H2O<real_t>() * M_H<real_t>() * (V*V) / m_H
            +
            // H2O*SO2 to HSO3 dissociation
            V * K_SO2<real_t>() / conc_H * SO2_helper
            +
            // HSO3 to SO3 dissociation 
            real_t(2) * // "2-" ion 
            V * K_SO2<real_t>() * K_HSO3<real_t>() / conc_H / conc_H * SO2_helper
            +
            // dissociation of S_VI to HSO4 (assumed there is no non-dissociated H2SO4)
            conc_H * m_S_VI / M_H2SO4<real_t>() / (conc_H + K_HSO4<real_t>())
            +
            // dissociation of HSO4 to SO4  (assumed there is no non-dissociated H2SO4)
            real_t(2) * // "2-" ion
            K_HSO4<real_t>() * m_S_VI / M_H2SO4<real_t>() / (conc_H + K_HSO4<real_t>())
            +
            // dissociation of CO2 * H2O to HCO3
            V * K_CO2<real_t>() / conc_H * CO2_helper
            +
            // dissociation of HCO3 to CO3
            real_t(2) * //"2-" ion
            V * K_CO2<real_t>() * K_HCO3<real_t>() / conc_H / conc_H * CO2_helper
            +
            // dissociation of HNO3 to NO3
            K_HNO3<real_t>() * n_HNO3_old / (K_HNO3<real_t>() + conc_H)
            - 
            // dissociation of NH3 * H2O to NH4
            K_NH3<real_t>() * n_NH3_old / (K_NH3<real_t>() + K_H2O<real_t>() / conc_H)
         )) / si::kilograms;
        } 
      };

      template <typename real_t>
      struct chem_electroneutral // TODO: does it have to be a struct/functor - perhaps ordinary function would suffice?
      { // uses toms748 scheme to solve for mass of H+ after dissociation
        // that satisfies electroneutrality
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t> &tpl) const
        {
          using namespace common::molar_mass;

          const quantity<si::amount, real_t>
            n_SO2_old  = thrust::get<0>(tpl) * si::moles,
            n_CO2_old  = thrust::get<1>(tpl) * si::moles,
            n_HNO3_old = thrust::get<2>(tpl) * si::moles,
            n_NH3_old  = thrust::get<3>(tpl) * si::moles;
          const quantity<si::mass, real_t>    m_S_VI  = thrust::get<4>(tpl) * si::kilograms;
          const quantity<si::volume, real_t>  V       = thrust::get<5>(tpl) * si::cubic_metres;
          
          // left side for search in toms748
          real_t m_H_pure = ((real_t(1e-7 * 1e3) * si::moles / si::cubic_metres) * V * M_H<real_t>()) / si::kilograms;

          uintmax_t max_iter = 44;

          real_t m_H = common::detail::toms748_solve(
	    detail::chem_minfun<real_t>(n_SO2_old, n_CO2_old, n_HNO3_old, n_NH3_old, m_S_VI, V),
            //m_H_pure, // min -> (pure water)
            real_t(1e-40),
	    real_t(1e-10), // max -> TODO
            common::detail::eps_tolerance<float>(sizeof(float) * 8), //TODO is it big enough?
            max_iter
	  ); 
          // std::cerr << "  " << real_t(1e-40) << " ... " << m_H << " ... " << "TODO" << real_t(1e-10) << std::endl;
          // TODO: asserts for K = f(m_H, m_...)
          return m_H;
        }
      };
 
      template <typename real_t, int chem_iter>
      struct chem_dissoc_diag // TODO: does it have to be a struct/functor - perhaps ordinary function would suffice?
      { // basing on the mass of H+ ions from chem_electroneutral
        // it recalculatess the mass of all other chemical species after dissociation
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::volume, real_t> 
            V      = thrust::get<0>(tpl) * si::cubic_metres;
          const quantity<si::mass, real_t>
            m_H    = thrust::get<1>(tpl) * si::kilograms,
            m_S_VI = thrust::get<2>(tpl) * si::kilograms;
          const quantity<si::amount, real_t>
            n_SO2_old  = thrust::get<3>(tpl) * si::moles,
            n_CO2_old  = thrust::get<4>(tpl) * si::moles,
            n_HNO3_old = thrust::get<5>(tpl) * si::moles,
            n_NH3_old  = thrust::get<6>(tpl) * si::moles;

          using namespace common::dissoc;     // K-prefixed
	  using namespace common::molar_mass; // M-prefixed

          quantity<common::amount_over_volume, real_t> conc_H, CO2_helper, SO2_helper;

          conc_H = m_H / M_H<real_t>() / V;
    
          CO2_helper = n_CO2_old / V 
                       / (real_t(1) + K_CO2<real_t>() / conc_H + K_CO2<real_t>() * K_HCO3<real_t>() / conc_H / conc_H);
          SO2_helper = n_SO2_old / V 
                       / (real_t(1) + K_SO2<real_t>() / conc_H + K_SO2<real_t>() * K_HSO3<real_t>() / conc_H / conc_H);

          switch (chem_iter)
          {
            case OH:
              return (
                M_H<real_t>() * M_OH<real_t>() 
                * K_H2O<real_t>()               // note: dissociation constant for pure water 
                * V * V / m_H                   // is actually k*[H2O] (Seinfeld and Pandis p 345)
              ) / si::kilograms; 
            case SO2:
              return (
                V * M_SO2_H2O<real_t>() * SO2_helper            
              ) / si::kilograms;
            case HSO3:
              return (
                V * M_HSO3<real_t>() * K_SO2<real_t>() / conc_H * SO2_helper
              ) / si::kilograms;
            case SO3:
              return (
                V * M_SO3<real_t>() * K_SO2<real_t>() * K_HSO3<real_t>() / conc_H / conc_H * SO2_helper
              ) / si::kilograms;
            case HSO4:
              return (
                 M_HSO4<real_t>() * m_S_VI / M_H2SO4<real_t>() * conc_H / (conc_H + K_HSO4<real_t>())
              ) / si::kilograms;
            case SO4:
              return (
                M_SO4<real_t>() * K_HSO4<real_t>() * m_S_VI / M_H2SO4<real_t>() / (conc_H + K_HSO4<real_t>())
              ) / si::kilograms;
            case CO2:
              return ( 
                V * M_CO2_H2O<real_t>() * CO2_helper
              ) / si::kilograms;
            case HCO3:
              return (
                V * M_HCO3<real_t>() * K_CO2<real_t>() / conc_H * CO2_helper
              ) / si::kilograms;
            case CO3:
              return ( 
                V * M_CO3<real_t>() * K_CO2<real_t>() * K_HCO3<real_t>() / conc_H / conc_H * CO2_helper
              ) / si::kilograms;
            case NH3:
              return (
                M_NH3_H2O<real_t>() * n_NH3_old * (K_H2O<real_t>() / conc_H) / (K_NH3<real_t>() + K_H2O<real_t>() / conc_H)
              ) / si::kilograms;
            case NH4:
              return (
                M_NH4<real_t>() * K_NH3<real_t>() * n_NH3_old / (K_NH3<real_t>() + K_H2O<real_t>() / conc_H)
              ) / si::kilograms;
            case HNO3:
              return ( 
                M_HNO3<real_t>() * n_HNO3_old * conc_H / (K_HNO3<real_t>() + conc_H) 
              ) / si::kilograms;
            case NO3:
              return (
                M_NO3<real_t>() * n_HNO3_old * K_HNO3<real_t>() / (K_HNO3<real_t>() + conc_H)
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
        assert(
          HNO3 == 0 && NH3  == 1 && CO2 == 2 &&
          SO2  == 3 && H2O2 == 4 && O3  == 5 
        );
        //Henry constant
        static const quantity<common::amount_over_volume_over_pressure, real_t> H_[chem_gas_n] = {
          H_HNO3<real_t>(), H_NH3<real_t>(),  H_CO2<real_t>(),
	  H_SO2<real_t>(),  H_H2O2<real_t>(), H_O3<real_t>()
        };
        //correction to Henry const due to temperature
        static const quantity<si::temperature, real_t> dHR_[chem_gas_n] = {
          dHR_HNO3<real_t>(), dHR_NH3<real_t>(),  dHR_CO2<real_t>(),
	  dHR_SO2<real_t>(),  dHR_H2O2<real_t>(), dHR_O3<real_t>()
        };
        //molar mass
        static const quantity<common::mass_over_amount, real_t> M_[chem_gas_n] = {
          M_HNO3<real_t>(),    M_NH3_H2O<real_t>(), M_CO2_H2O<real_t>(),
          M_SO2_H2O<real_t>(), M_H2O2<real_t>(), M_O3<real_t>()
        };
        //gas phase diffusivity
        static const quantity<common::diffusivity, real_t> D_[chem_gas_n] = {
          D_HNO3<real_t>(), D_NH3<real_t>(),  D_CO2<real_t>(),
          D_SO2<real_t>(),  D_H2O2<real_t>(), D_O3<real_t>()
        };
        //accomodation coefficient
         static const quantity<si::dimensionless, real_t> ac_[chem_gas_n] = {
          ac_HNO3<real_t>(), ac_NH3<real_t>(),  ac_CO2<real_t>(),
          ac_SO2<real_t>(),  ac_H2O2<real_t>(), ac_O3<real_t>()
        };

	typedef thrust::permutation_iterator<
	  typename thrust_device::vector<real_t>::iterator,
	  typename thrust_device::vector<thrust_size_t>::iterator
	> pi_t;

	typedef thrust::zip_iterator<
          thrust::tuple<
            pi_t,                                              // pressure
            pi_t,                                              // temprarature
            typename thrust_device::vector<real_t>::iterator,  // m_old
            typename thrust_device::vector<real_t>::iterator   // rw2
          >
        > zip_it_t;

        for (int i = 0; i < chem_gas_n; ++i)
        {
          thrust::transform(
            V.begin(), V.end(),            // input - 1st arg
            zip_it_t(thrust::make_tuple(   // input - 2nd arg
	      thrust::make_permutation_iterator(p.begin(), ijk.begin()),
	      thrust::make_permutation_iterator(T.begin(), ijk.begin()),
              chem_bgn[i],
              rw2.begin() 
            )),
	    chem_bgn[i],                                                                                        // output
	    detail::chem_Henry_fun<real_t>(H_[i], dHR_[i], M_[i], chem_gas[i], D_[i], ac_[i], dt * si::seconds) // op
	  );
        }
      }

      if (chem_dsc == true){  //TODO move to a separate function and then move the logic to opts particle_step 
        // 2/4: equilibrium stuff: dissociation
 
        // save number of moles of dissolving chem species
        // could be avoided if we could have more than 10 elements in a tuple...
        // TODO - add thrust tuple for 12 arguments https://github.com/thrust/thrust
        thrust_device::vector<real_t> &n_SO2_old(tmp_device_real_part_SO2);
        thrust_device::vector<real_t> &n_CO2_old(tmp_device_real_part_CO2);
        thrust_device::vector<real_t> &n_HNO3_old(tmp_device_real_part_HNO3);
        thrust_device::vector<real_t> &n_NH3_old(tmp_device_real_part_NH3);

        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::const_iterator, // SO2   CO2
            typename thrust_device::vector<real_t>::const_iterator, // HSO3  HCO3
            typename thrust_device::vector<real_t>::const_iterator  // SO3   CO3
          >
        > zip_it_t_3;

        thrust::transform(
          zip_it_t_3(thrust::make_tuple(chem_bgn[SO2], chem_bgn[HSO3], chem_bgn[SO3])),  //input begin
          zip_it_t_3(thrust::make_tuple(chem_end[SO2], chem_end[HSO3], chem_end[SO3])),  //input end
          n_SO2_old.begin(),                                                          //output
          detail::chem_save_moles_2<real_t>(M_SO2_H2O<real_t>(), M_HSO3<real_t>(), M_SO3<real_t>()) //op
        );

        thrust::transform(
          zip_it_t_3(thrust::make_tuple(chem_bgn[CO2], chem_bgn[HCO3], chem_bgn[CO3])),  //input begin
          zip_it_t_3(thrust::make_tuple(chem_end[CO2], chem_end[HCO3], chem_end[CO3])),  //input end
          n_CO2_old.begin(),                                                          //output
          detail::chem_save_moles_2<real_t>(M_CO2_H2O<real_t>(), M_HCO3<real_t>(), M_CO3<real_t>()) //op
        );

        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::const_iterator, // HNO3  NH3
            typename thrust_device::vector<real_t>::const_iterator  // NO3   NH4
          >
        > zip_it_t_2;

        thrust::transform(
          zip_it_t_2(thrust::make_tuple(chem_bgn[HNO3], chem_bgn[NO3])),  //input begin
          zip_it_t_2(thrust::make_tuple(chem_end[HNO3], chem_end[NO3])),  //input end
          n_HNO3_old.begin(),                                           //output
          detail::chem_save_moles_1<real_t>(M_HNO3<real_t>(), M_NO3<real_t>()) //op
        );

        thrust::transform(
          zip_it_t_2(thrust::make_tuple(chem_bgn[NH3], chem_bgn[NH4])),  //input begin
          zip_it_t_2(thrust::make_tuple(chem_end[NH3], chem_end[NH4])),  //input end
          n_NH3_old.begin(),                                           //output
          detail::chem_save_moles_1<real_t>(M_NH3_H2O<real_t>(), M_NH4<real_t>()) //op
        );

        { // calculate H+ ions after dissociation so that drops remain electroneutral
          typedef thrust::zip_iterator<
            thrust::tuple<
              typename thrust_device::vector<real_t>::iterator, // SO2_old
              typename thrust_device::vector<real_t>::iterator, // CO2_old
              typename thrust_device::vector<real_t>::iterator, // HNO3_old
              typename thrust_device::vector<real_t>::iterator, // NH3_old
              typename thrust_device::vector<real_t>::iterator, // S_VI_old
              typename thrust_device::vector<real_t>::iterator  // V
            >
          > zip_it_t;

          thrust::transform(
	    zip_it_t(thrust::make_tuple(
              n_SO2_old.begin(),n_CO2_old.begin(), n_HNO3_old.begin(), n_NH3_old.begin(),
              chem_bgn[S_VI], 
              V.begin())),  // input - begin
	    zip_it_t(thrust::make_tuple(
              n_SO2_old.end(), n_CO2_old.end(), n_HNO3_old.end(), n_NH3_old.end(),
              chem_end[S_VI], 
              V.end())),  // input - end
	    chem_bgn[H],                                                             // output
	    detail::chem_electroneutral<real_t>()                                    // op
	  );
        }

        { // diagnose the rest of ions basing on H+
          typedef thrust::zip_iterator<
            thrust::tuple<
              typename thrust_device::vector<real_t>::iterator, // V
              typename thrust_device::vector<real_t>::iterator, // H 
              typename thrust_device::vector<real_t>::iterator, // m_S_VI
              typename thrust_device::vector<real_t>::iterator, // n_SO2_old
              typename thrust_device::vector<real_t>::iterator, // n_CO2_old
              typename thrust_device::vector<real_t>::iterator, // n_HNO3_old
              typename thrust_device::vector<real_t>::iterator  // n_NH3_old
            >
          > zip_it_t;

          zip_it_t 
            arg_begin(thrust::make_tuple(V.begin(), 
              chem_bgn[H], chem_bgn[S_VI], n_SO2_old.begin(), n_CO2_old.begin(), n_HNO3_old.begin(), n_NH3_old.begin())
            ),
            arg_end(thrust::make_tuple(V.end(),
              chem_end[H], chem_end[S_VI], n_SO2_old.begin(), n_CO2_old.begin(), n_HNO3_old.begin(), n_NH3_old.begin())
            );
        
          thrust::transform(arg_begin, arg_end, chem_bgn[OH],   detail::chem_dissoc_diag<real_t, OH  >());
          thrust::transform(arg_begin, arg_end, chem_bgn[SO2],  detail::chem_dissoc_diag<real_t, SO2 >()); 
          thrust::transform(arg_begin, arg_end, chem_bgn[HSO3], detail::chem_dissoc_diag<real_t, HSO3>()); 
          thrust::transform(arg_begin, arg_end, chem_bgn[SO3],  detail::chem_dissoc_diag<real_t, SO3 >());
          thrust::transform(arg_begin, arg_end, chem_bgn[HSO4], detail::chem_dissoc_diag<real_t, HSO4>());
          thrust::transform(arg_begin, arg_end, chem_bgn[SO4],  detail::chem_dissoc_diag<real_t, SO4 >());
          thrust::transform(arg_begin, arg_end, chem_bgn[CO2],  detail::chem_dissoc_diag<real_t, CO2 >()); 
          thrust::transform(arg_begin, arg_end, chem_bgn[HCO3], detail::chem_dissoc_diag<real_t, HCO3>()); 
          thrust::transform(arg_begin, arg_end, chem_bgn[CO3],  detail::chem_dissoc_diag<real_t, CO3 >());
          thrust::transform(arg_begin, arg_end, chem_bgn[HNO3], detail::chem_dissoc_diag<real_t, HNO3>());
          thrust::transform(arg_begin, arg_end, chem_bgn[NO3],  detail::chem_dissoc_diag<real_t, NO3 >());
          thrust::transform(arg_begin, arg_end, chem_bgn[NH3],  detail::chem_dissoc_diag<real_t, NH3 >());
          thrust::transform(arg_begin, arg_end, chem_bgn[NH4],  detail::chem_dissoc_diag<real_t, NH4 >());
/*
          std::cerr<<" "<<std::endl;
          std::cerr<<"positive ions: "<< std::endl;
          std::cerr<<"H+"<<std::endl;
          debug::print(chem_bgn[H], chem_end[H]);
          std::cerr<<"NH4+"<<std::endl;
          debug::print(chem_bgn[NH4], chem_end[NH4]);
          std::cerr<<"negative ions: "<< std::endl;

          std::cerr<<"OH-"<<std::endl;
          debug::print(chem_bgn[OH], chem_end[OH]);

          std::cerr<<"HSO3-"<<std::endl;
          debug::print(chem_bgn[HSO3], chem_end[HSO3]);
          std::cerr<<"SO3--"<<std::endl;
          debug::print(chem_bgn[SO3], chem_end[SO3]);

          std::cerr<<"HSO4-"<<std::endl;
          debug::print(chem_bgn[HSO4], chem_end[HSO4]);
          std::cerr<<"SO4--"<<std::endl;
          debug::print(chem_bgn[SO4], chem_end[SO4]);

          std::cerr<<"HCO3-"<<std::endl;
          debug::print(chem_bgn[HCO3], chem_end[HCO3]);
          std::cerr<<"CO3--"<<std::endl;
          debug::print(chem_bgn[CO3], chem_end[CO3]);

          std::cerr<<"NO3-"<<std::endl;
          debug::print(chem_bgn[NO3], chem_end[NO3]);
*/        }
      }

      if (chem_rct == true){  //TODO move to a separate function and then move the logic to opts particle_step 
        // 3/4: non-equilibrium stuff
        thrust_device::vector<real_t> &old_S_VI(tmp_device_real_part_2);
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

        // 4/4: recomputing dry radii
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
    }
  };  
};
