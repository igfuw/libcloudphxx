// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <libcloudph++/common/detail/toms748.hpp>
#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/dissoc.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
/*
      template <typename real_t>
      struct chem_save_moles_2
      {
        // helper for dissociation - saves the number of moles of dissolving chem species
        // could be avoided if we could have more than 10 elements in a tuple...
        // TODO - add thrust tuple for 12 arguments https://github.com/thrust/thrust
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
        // helper for dissociation - saves the number of moles of dissolving chem species
        // could be avoided if we could have more than 10 elements in a tuple...
        // TODO - add thrust tuple for 12 arguments https://github.com/thrust/thrust
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
 */
      template <typename real_t>
      struct chem_minfun
      { // helper function for maintaining electroneutrality
        // returns the difference between the mass of H+ ions already present 
        // and the mass of H+ ions needed to be added (due to dissociation)
        // in order to maintain electroneutrality
	const quantity<si::mass, real_t> m_S_IV, m_C_IV, m_N_III, m_N_V, m_S_VI;
	const quantity<si::volume, real_t> V;
        const quantity<si::temperature, real_t> T; 
        
        // ctor
        BOOST_GPU_ENABLED
        chem_minfun(
	  const quantity<si::mass, real_t> &m_S_IV,
	  const quantity<si::mass, real_t> &m_C_IV,
	  const quantity<si::mass, real_t> &m_N_V,
	  const quantity<si::mass, real_t> &m_N_III,
	  const quantity<si::mass, real_t> &m_S_VI,
	  const quantity<si::volume, real_t> &V,
          const quantity<si::temperature, real_t> T
        ) :
          m_S_IV(m_S_IV), m_C_IV(m_C_IV), m_N_V(m_N_V), m_N_III(m_N_III), 
          m_S_VI(m_S_VI), V(V), T(T)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &arg) const
        {
	  using namespace common::molar_mass;
	  using namespace common::dissoc;
          
          //calculate temperature dependant dissociation constants
          // TODO repeated in chem dissoc diag
          quantity<common::amount_over_volume, real_t> Kt_CO2, Kt_HCO3, Kt_SO2, Kt_HSO3, Kt_NH3, Kt_HNO3, Kt_HSO4;
          Kt_CO2  = K_temp(T, K_CO2<real_t>(),  dKR_CO2<real_t>());
          Kt_HCO3 = K_temp(T, K_HCO3<real_t>(), dKR_HCO3<real_t>());
          Kt_SO2  = K_temp(T, K_SO2<real_t>(),  dKR_SO2<real_t>());
          Kt_HSO3 = K_temp(T, K_HSO3<real_t>(), dKR_HSO3<real_t>());
          Kt_NH3  = K_temp(T, K_NH3<real_t>(),  dKR_NH3<real_t>());
          Kt_HNO3 = K_temp(T, K_HNO3<real_t>(), dKR_HNO3<real_t>());
          Kt_HSO4 = K_temp(T, K_HSO4<real_t>(), dKR_HSO4<real_t>());

          //helper for concentration of H+ ions
          const quantity<si::mass, real_t> m_H = arg * si::kilograms;
          quantity<common::amount_over_volume, real_t> conc_H;
          conc_H = m_H / M_H<real_t>() / V;

          return (-m_H + M_H<real_t>() * (
            // dissociation of pure water 
            K_H2O<real_t>() * M_H<real_t>() * (V*V) / m_H
            +
            // H2O*SO2 to HSO3 dissociation
            m_S_IV / M_SO2_H2O<real_t>() * Kt_SO2 / conc_H 
              / (real_t(1) + Kt_SO2 / conc_H + Kt_SO2 * Kt_HSO3 / conc_H / conc_H)
            +
            // HSO3 to SO3 dissociation 
            real_t(2) * // "2-" ion 
            m_S_IV / M_SO2_H2O<real_t>() * Kt_SO2 * Kt_HSO3 / conc_H / conc_H 
              / (real_t(1) + Kt_SO2 / conc_H + Kt_SO2 * Kt_HSO3 / conc_H / conc_H)
            +
            // dissociation of S_VI to HSO4 (assumed there is no non-dissociated H2SO4)
            conc_H * m_S_VI / M_H2SO4<real_t>() / (conc_H + Kt_HSO4)
            +
            // dissociation of HSO4 to SO4  (assumed there is no non-dissociated H2SO4)
            real_t(2) * // "2-" ion
            Kt_HSO4 * m_S_VI / M_H2SO4<real_t>() / (conc_H + Kt_HSO4)
            +
            // dissociation of CO2 * H2O to HCO3
            m_C_IV / M_CO2_H2O<real_t>() * Kt_CO2 / conc_H 
              / (real_t(1) + Kt_CO2 / conc_H + Kt_CO2 * Kt_HCO3 / conc_H / conc_H)
            +
            // dissociation of HCO3 to CO3
            real_t(2) * //"2-" ion
            m_C_IV / M_CO2_H2O<real_t>() * Kt_CO2 * Kt_HCO3 / conc_H / conc_H 
              / (real_t(1) + Kt_CO2 / conc_H + Kt_CO2 * Kt_HCO3 / conc_H / conc_H)
            +
            // dissociation of HNO3 to NO3
            m_N_V / M_HNO3<real_t>() * Kt_HNO3 / conc_H / (real_t(1.) + Kt_HNO3 / conc_H)
            - 
            // dissociation of NH3 * H2O to NH4
            m_N_III / M_NH3_H2O<real_t>() * Kt_NH3 / K_H2O<real_t>() * conc_H / (real_t(1) + Kt_NH3 / K_H2O<real_t>() * conc_H)
         )) / si::kilograms;
        } 
      };

      template <typename real_t>
      struct chem_electroneutral // TODO: does it have to be a struct/functor - perhaps ordinary function would suffice?
      { // uses toms748 scheme to solve for mass of H+ after dissociation
        // that satisfies electroneutrality
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t> &tpl) const
        {
          using namespace common::molar_mass;

          const quantity<si::mass, real_t>
            m_S_IV  = thrust::get<0>(tpl) * si::kilograms,
            m_C_IV  = thrust::get<1>(tpl) * si::kilograms,
            m_N_V   = thrust::get<2>(tpl) * si::kilograms,
            m_N_III = thrust::get<3>(tpl) * si::kilograms,
            m_S_VI  = thrust::get<4>(tpl) * si::kilograms;
          const quantity<si::volume, real_t>  V       = thrust::get<5>(tpl) * si::cubic_metres;
          const quantity<si::temperature, real_t> T   = thrust::get<6>(tpl) * si::kelvins;
          
          // limits for search in toms748
          real_t m_H_rht = ((real_t(1e1  * 1e3)  * si::moles / si::cubic_metres) * V * M_H<real_t>()) / si::kilograms;
          real_t m_H_lft = ((real_t(1e-9 * 1e3) * si::moles / si::cubic_metres) * V * M_H<real_t>()) / si::kilograms;

          uintmax_t max_iter = 44;

          real_t m_H = common::detail::toms748_solve(
	    detail::chem_minfun<real_t>(m_S_IV, m_C_IV, m_N_V, m_N_III, m_S_VI, V, T),
            m_H_lft,
	    m_H_rht,
            common::detail::eps_tolerance<float>(sizeof(float) * 8), //TODO is it big enough?
            max_iter
	  ); 

          //real_t ph_helper = real_t(-1.) * log10(m_H / (M_H<real_t>() / si::kilograms * si::moles) / (V / si::cubic_metres) / real_t(1000.));
          //std::cerr << "  " << m_H_lft << " ... " << m_H << " ... " << m_H_rht << " -> "<< ph_helper<< std::endl;
          // TODO: asserts for K = f(m_H, m_...)
          return m_H;
        }
      };
 
/*
      template <typename real_t, int chem_iter>
      struct chem_dissoc_diag // TODO: does it have to be a struct/functor - perhaps ordinary function would suffice?
      { // basing on the mass of H+ ions from chem_electroneutral
        // it recalculatess the mass of all other chemical species after dissociation
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t> &tpl) const
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
          const quantity<si::temperature, real_t>
            T      = thrust::get<7>(tpl) * si::kelvins;

          using namespace common::dissoc;     // K-prefixed
	  using namespace common::molar_mass; // M-prefixed

          //calculate temperature dependant dissociation constants
          // TODO - repeated in chem_min_fun
          quantity<common::amount_over_volume, real_t> Kt_CO2, Kt_HCO3, Kt_SO2, Kt_HSO3, Kt_NH3, Kt_HNO3, Kt_HSO4;
          Kt_CO2  = K_temp(T, K_CO2<real_t>(),  dKR_CO2<real_t>());
          Kt_HCO3 = K_temp(T, K_HCO3<real_t>(), dKR_HCO3<real_t>());
          Kt_SO2  = K_temp(T, K_SO2<real_t>(),  dKR_SO2<real_t>());
          Kt_HSO3 = K_temp(T, K_HSO3<real_t>(), dKR_HSO3<real_t>());
          Kt_NH3  = K_temp(T, K_NH3<real_t>(),  dKR_NH3<real_t>());
          Kt_HNO3 = K_temp(T, K_HNO3<real_t>(), dKR_HNO3<real_t>());
          Kt_HSO4 = K_temp(T, K_HSO4<real_t>(),  dKR_HSO4<real_t>());

          quantity<common::amount_over_volume, real_t> conc_H, CO2_helper, SO2_helper;

          conc_H = m_H / M_H<real_t>() / V;
    
          CO2_helper = n_CO2_old / V 
                       / (real_t(1) + Kt_CO2 / conc_H + Kt_CO2 * Kt_HCO3 / conc_H / conc_H);
          SO2_helper = n_SO2_old / V 
                       / (real_t(1) + Kt_SO2 / conc_H + Kt_SO2 * Kt_HSO3 / conc_H / conc_H);

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
                V * M_HSO3<real_t>() * Kt_SO2 / conc_H * SO2_helper
              ) / si::kilograms;
            case SO3:
              return (
                V * M_SO3<real_t>() * Kt_SO2 * Kt_HSO3 / conc_H / conc_H * SO2_helper
              ) / si::kilograms;
            case HSO4:
              return (
                 M_HSO4<real_t>() * m_S_VI / M_H2SO4<real_t>() * conc_H / (conc_H + Kt_HSO4)
              ) / si::kilograms;
            case SO4:
              return (
                M_SO4<real_t>() * Kt_HSO4 * m_S_VI / M_H2SO4<real_t>() / (conc_H + Kt_HSO4)
              ) / si::kilograms;
            case CO2:
              return ( 
                V * M_CO2_H2O<real_t>() * CO2_helper
              ) / si::kilograms;
            case HCO3:
              return (
                V * M_HCO3<real_t>() * Kt_CO2 / conc_H * CO2_helper
              ) / si::kilograms;
            case CO3:
              return ( 
                V * M_CO3<real_t>() * Kt_CO2 * Kt_HCO3 / conc_H / conc_H * CO2_helper
              ) / si::kilograms;
            case NH3:
              return (
                M_NH3_H2O<real_t>() * n_NH3_old * (K_H2O<real_t>() / conc_H) / (Kt_NH3 + K_H2O<real_t>() / conc_H)
              ) / si::kilograms;
            case NH4:
              return (
                M_NH4<real_t>() * Kt_NH3 * n_NH3_old / (Kt_NH3 + K_H2O<real_t>() / conc_H)
              ) / si::kilograms;
            case HNO3:
              return ( 
                M_HNO3<real_t>() * n_HNO3_old * conc_H / (Kt_HNO3 + conc_H) 
              ) / si::kilograms;
            case NO3:
              return (
                M_NO3<real_t>() * n_HNO3_old * Kt_HNO3 / (Kt_HNO3 + conc_H)
              ) / si::kilograms;
            default:
              assert(false);
              return 0;
          }
        }
      };
    };
*/
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_dissoc()
    {   
      using namespace common::molar_mass; // M-prefixed

      thrust_device::vector<real_t> &V(tmp_device_real_part);

      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

/* 
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
        n_SO2_old.begin(),                                                             //output
        detail::chem_save_moles_2<real_t>(M_SO2_H2O<real_t>(), M_HSO3<real_t>(), M_SO3<real_t>()) //op
      );

      thrust::transform(
        zip_it_t_3(thrust::make_tuple(chem_bgn[CO2], chem_bgn[HCO3], chem_bgn[CO3])),  //input begin
        zip_it_t_3(thrust::make_tuple(chem_end[CO2], chem_end[HCO3], chem_end[CO3])),  //input end
        n_CO2_old.begin(),                                                             //output
        detail::chem_save_moles_2<real_t>(M_CO2_H2O<real_t>(), M_HCO3<real_t>(), M_CO3<real_t>()) //op
      );

      typedef thrust::zip_iterator<
        thrust::tuple<
          typename thrust_device::vector<real_t>::const_iterator, // HNO3  NH3
          typename thrust_device::vector<real_t>::const_iterator  // NO3   NH4
        >
      > zip_it_t_2;

      thrust::transform(
        zip_it_t_2(thrust::make_tuple(chem_bgn[HNO3], chem_bgn[NO3])),       //input begin
        zip_it_t_2(thrust::make_tuple(chem_end[HNO3], chem_end[NO3])),       //input end
        n_HNO3_old.begin(),                                                  //output
        detail::chem_save_moles_1<real_t>(M_HNO3<real_t>(), M_NO3<real_t>()) //op
      );

      thrust::transform(
        zip_it_t_2(thrust::make_tuple(chem_bgn[NH3], chem_bgn[NH4])),           //input begin
        zip_it_t_2(thrust::make_tuple(chem_end[NH3], chem_end[NH4])),           //input end
        n_NH3_old.begin(),                                                      //output
        detail::chem_save_moles_1<real_t>(M_NH3_H2O<real_t>(), M_NH4<real_t>()) //op
      );

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_t;
*/
      { // calculate H+ ions after dissociation so that drops remain electroneutral
        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::iterator, // S_IV
            typename thrust_device::vector<real_t>::iterator, // C_IV
            typename thrust_device::vector<real_t>::iterator, // N_V
            typename thrust_device::vector<real_t>::iterator, // N_III
            typename thrust_device::vector<real_t>::iterator, // S_VI
            typename thrust_device::vector<real_t>::iterator, // V
            pi_t  // T
          >
        > zip_it_t;

        thrust::transform(
          zip_it_t(thrust::make_tuple(
            chem_bgn[SO2], chem_bgn[CO2], chem_bgn[HNO3], chem_bgn[NH3], chem_bgn[S_VI], 
            V.begin(),
            thrust::make_permutation_iterator(T.begin(), ijk.begin()))),                   // input - begin
          zip_it_t(thrust::make_tuple(
            chem_end[SO2], chem_end[CO2], chem_end[HNO3], chem_end[NH3], chem_end[S_VI], 
            V.end(),
            thrust::make_permutation_iterator(T.end(), ijk.end()))),                       // input - end
         chem_bgn[H],                                                                      // output
         detail::chem_electroneutral<real_t>()                                             // op
        );
      }
/*
      { // diagnose the rest of ions basing on H+
        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::iterator, // V
            typename thrust_device::vector<real_t>::iterator, // H 
            typename thrust_device::vector<real_t>::iterator, // m_S_VI
            typename thrust_device::vector<real_t>::iterator, // n_SO2_old
            typename thrust_device::vector<real_t>::iterator, // n_CO2_old
            typename thrust_device::vector<real_t>::iterator, // n_HNO3_old
            typename thrust_device::vector<real_t>::iterator, // n_NH3_old
            pi_t  // T
          >
        > zip_it_t;

        zip_it_t 
          arg_begin(thrust::make_tuple(V.begin(), 
            chem_bgn[H], chem_bgn[S_VI], n_SO2_old.begin(), 
            n_CO2_old.begin(), n_HNO3_old.begin(), n_NH3_old.begin(), 
            thrust::make_permutation_iterator(T.begin(), ijk.begin())
        )),
          arg_end(thrust::make_tuple(V.end(),
            chem_end[H], chem_end[S_VI], n_SO2_old.end(),
            n_CO2_old.end(), n_HNO3_old.end(), n_NH3_old.end(),
            thrust::make_permutation_iterator(T.end(), ijk.end())
        ));
        
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
      }
*/
      for (int i = 0; i < chem_gas_n; ++i){
        assert(boost::math::isfinite(*thrust::min_element(chem_bgn[i], chem_end[i])));
      }
    }
  };  
};
