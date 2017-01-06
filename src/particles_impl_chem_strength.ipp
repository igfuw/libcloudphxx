// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <class real_t>
      struct set_chem_flag
      { // helper for choosing which droplets are dilute and can have chemical reactions
        // calculates the ionic strength of droplets
        BOOST_GPU_ENABLED        
        unsigned int operator()(const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t> &tpl) const 
        {
          using namespace common::molar_mass;
          using namespace common::dissoc;

          //handle tuple elements
          const quantity<si::volume, real_t>  V       = thrust::get<6>(tpl) * si::cubic_metres;
          const quantity<si::temperature, real_t> T   = thrust::get<7>(tpl) * si::kelvins;
          const quantity<si::area, real_t> rw2        = thrust::get<8>(tpl) * si::square_metres;
          const quantity<common::amount_over_volume, real_t>
            conc_S_IV  = thrust::get<0>(tpl) * si::kilograms / M_SO2_H2O<real_t>() / V;
          const quantity<common::amount_over_volume, real_t>
            conc_C_IV  = thrust::get<1>(tpl) * si::kilograms / M_CO2_H2O<real_t>() / V;
          const quantity<common::amount_over_volume, real_t>
            conc_N_V   = thrust::get<2>(tpl) * si::kilograms / M_HNO3<real_t>()    / V;
          const quantity<common::amount_over_volume, real_t>
            conc_N_III = thrust::get<3>(tpl) * si::kilograms / M_NH3_H2O<real_t>() / V;
          const quantity<common::amount_over_volume, real_t>
            conc_S_VI  = thrust::get<4>(tpl) * si::kilograms / M_H2SO4<real_t>()   / V;
          const quantity<common::amount_over_volume, real_t>
            conc_H     = thrust::get<5>(tpl) * si::kilograms / M_H<real_t>()       / V;

          //calculate temperature dependant dissociation constants - TODO - copied from dissoc!
          quantity<common::amount_over_volume, real_t> Kt_CO2  = K_temp(T, K_CO2<real_t>(),  dKR_CO2<real_t>());
          quantity<common::amount_over_volume, real_t> Kt_HCO3 = K_temp(T, K_HCO3<real_t>(), dKR_HCO3<real_t>());
          quantity<common::amount_over_volume, real_t> Kt_SO2  = K_temp(T, K_SO2<real_t>(),  dKR_SO2<real_t>());
          quantity<common::amount_over_volume, real_t> Kt_HSO3 = K_temp(T, K_HSO3<real_t>(), dKR_HSO3<real_t>());
          quantity<common::amount_over_volume, real_t> Kt_NH3  = K_temp(T, K_NH3<real_t>(),  dKR_NH3<real_t>());
          quantity<common::amount_over_volume, real_t> Kt_HNO3 = K_temp(T, K_HNO3<real_t>(), dKR_HNO3<real_t>());
          quantity<common::amount_over_volume, real_t> Kt_HSO4 = K_temp(T, K_HSO4<real_t>(), dKR_HSO4<real_t>());
 
          return 
            //(rw2 / si::square_metres > 4e-12 
            //  && 
              ((real_t(0.5) * (
                conc_H * conc_S_VI / (conc_H + Kt_HSO4) +
                real_t(4) * Kt_HSO4 * conc_S_VI / (conc_H + Kt_HSO4) +
                Kt_CO2 * conc_H  * conc_C_IV / (conc_H * conc_H + Kt_CO2 * conc_H + Kt_CO2 * Kt_HCO3) + 
                real_t(4) * Kt_CO2 * Kt_HCO3 * conc_C_IV / (conc_H * conc_H + Kt_CO2 * conc_H + Kt_CO2 * Kt_HCO3) +
                Kt_SO2 * conc_H * conc_S_IV / (conc_H * conc_H + Kt_SO2 * conc_H + Kt_SO2 * Kt_HSO3) +
                real_t(4) * Kt_SO2 * Kt_HSO3 * conc_S_IV / (conc_H * conc_H + Kt_SO2 * conc_H + Kt_SO2 * Kt_HSO3) +
                Kt_HNO3 * conc_N_V / (conc_H + Kt_HNO3) + 
                Kt_NH3 * conc_H * conc_N_III / (K_H2O<real_t>() + Kt_NH3 * conc_H)
              ) * si::cubic_metres / si::moles) < 0.02 * 1000)//) 
            ? 1 : 0;                                  // ^^^^ convert from moles/litr to moles/m3
        }                                            
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_flag_ante()
    { 
      thrust_device::vector<unsigned int> &chem_flag(tmp_device_n_part);
      thrust_device::vector<real_t> &V(tmp_device_real_part);

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_t;

      {
        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::iterator, // S_IV
            typename thrust_device::vector<real_t>::iterator, // C_IV
            typename thrust_device::vector<real_t>::iterator, // N_V
            typename thrust_device::vector<real_t>::iterator, // N_III
            typename thrust_device::vector<real_t>::iterator, // S_VI
            typename thrust_device::vector<real_t>::iterator, // H
            typename thrust_device::vector<real_t>::iterator, // V
            pi_t,                                             // T
            typename thrust_device::vector<real_t>::iterator  // rw2
          >
        > zip_it_t;

        // set flag to those SDs that take part in chemical reactions
        thrust::transform(
          zip_it_t(thrust::make_tuple(
            chem_bgn[SO2], chem_bgn[CO2], chem_bgn[HNO3], chem_bgn[NH3], chem_bgn[S_VI], chem_bgn[H],
            V.begin(), thrust::make_permutation_iterator(T.begin(), ijk.begin()), rw2.begin())),                   // input - begin
          zip_it_t(thrust::make_tuple(
            chem_end[SO2], chem_end[CO2], chem_end[HNO3], chem_end[NH3], chem_end[S_VI], chem_end[H],
            V.end(), thrust::make_permutation_iterator(T.end(), ijk.end()), rw2.end())),                           // input - end
          chem_flag.begin(),               // output
          detail::set_chem_flag<real_t>()  // op
        );
      }
    }
  };  
};
