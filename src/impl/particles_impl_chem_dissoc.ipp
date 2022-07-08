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
          real_t m_H_rht = ((real_t(1e1  * 1e3) * si::moles / si::cubic_metres) * V * M_H<real_t>()) / si::kilograms;
          real_t m_H_lft = ((real_t(1e-8 * 1e3) * si::moles / si::cubic_metres) * V * M_H<real_t>()) / si::kilograms;

          uintmax_t max_iter = 100;

          real_t m_H = common::detail::toms748_solve(
	    detail::chem_minfun<real_t>(m_S_IV, m_C_IV, m_N_V, m_N_III, m_S_VI, V, T),
            m_H_lft,
	    m_H_rht,
            common::detail::eps_tolerance<float>(sizeof(float) * 8), //TODO is it big enough?
            max_iter
	  ); 
/*
          real_t ph_helper = real_t(-1.) * log10(m_H / (M_H<real_t>() / si::kilograms * si::moles) 
                             / (V / si::cubic_metres) / real_t(1000.));
          std::cerr << "  " << m_H_lft << " ... " << m_H << " ... " << m_H_rht << " -> pH  = "<< ph_helper<< std::endl;
          // TODO: asserts for K = f(m_H, m_...)
*/
          return m_H;
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_dissoc()
    {   
      using namespace common::molar_mass; // M-prefixed

      thrust_device::vector<real_t> &V(tmp_device_real_part);
      const thrust_device::vector<unsigned int> &chem_flag(tmp_device_n_part);

      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_t;

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

        thrust::transform_if(
          zip_it_t(thrust::make_tuple(
            chem_bgn[SO2], chem_bgn[CO2], chem_bgn[HNO3], chem_bgn[NH3], chem_bgn[S_VI], 
            V.begin(),
            thrust::make_permutation_iterator(T_ref.begin(), ijk.begin_ref()))),                   // input - begin
          zip_it_t(thrust::make_tuple(
            chem_end[SO2], chem_end[CO2], chem_end[HNO3], chem_end[NH3], chem_end[S_VI], 
            V.end(),
            thrust::make_permutation_iterator(T_ref.end(), ijk.end_ref()))),                       // input - end
         chem_flag.begin(),                                                                // stencil
         chem_bgn[H],                                                                      // output
         detail::chem_electroneutral<real_t>(),                                            // op
         thrust::identity<unsigned int>()
        );
      }

#if !defined(__NVCC__)
      using boost::math::isfinite;
#endif
      for (int i = 0; i < chem_gas_n; ++i){
        assert(isfinite(*thrust::min_element(chem_bgn[i], chem_end[i])));
      }
    }
  };  
};
