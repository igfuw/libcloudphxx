// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/henry.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename n_t, typename real_t>
      struct chem_summator
      { // calculate the mass of chem compounds (multiplicity * mass)  
        BOOST_GPU_ENABLED
        real_t operator()(const n_t &n, const real_t &chem) const
        {
          return n * chem;
        }
      };

      template <typename real_t>
      struct chem_Henry_fun
      { // gas absorption into cloud droplets (Henrys law)
        const quantity<common::amount_over_volume_over_pressure, real_t> H;
        const quantity<si::temperature, real_t> dHR;
        const quantity<common::mass_over_amount, real_t> M_gas;
        const quantity<common::mass_over_amount, real_t> M_aq;
        const quantity<common::diffusivity, real_t> D;
        const quantity<si::dimensionless, real_t> acc_coeff;
        const quantity<si::time, real_t> dt;

        // ctor
        chem_Henry_fun(
          const quantity<common::amount_over_volume_over_pressure, real_t> &H,
          const quantity<si::temperature, real_t> &dHR,
          const quantity<common::mass_over_amount, real_t> &M_gas,
          const quantity<common::mass_over_amount, real_t> &M_aq,
          const quantity<common::diffusivity, real_t> &D,
          const quantity<si::dimensionless, real_t> &acc_coeff,
          const quantity<si::time, real_t> &dt
        ) : 
          H(H), dHR(dHR), M_gas(M_gas), M_aq(M_aq), D(D), acc_coeff(acc_coeff), dt(dt)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(
          const real_t &V, 
          const thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::pressure, real_t>      p      = thrust::get<0>(tpl) * si::pascals; 
          const quantity<si::temperature, real_t>   T      = thrust::get<1>(tpl) * si::kelvins;     
          const quantity<si::dimensionless, real_t> c      = thrust::get<2>(tpl);     
          const quantity<si::mass, real_t>          m_old  = thrust::get<3>(tpl) * si::kilograms;     
          const quantity<si::area, real_t>          rw2    = thrust::get<4>(tpl) * si::metres * si::metres; 
          const quantity<si::mass_density, real_t>  rhod   = thrust::get<5>(tpl) * si::kilograms / si::cubic_metres; 
          const quantity<si::volume, real_t>        V_old  = thrust::get<6>(tpl) * si::cubic_metres;

          //implicit solution to the eq. 8.22 from chapter 8.4.2 in Peter Warneck Chemistry of the Natural Atmosphere  
          return  ((m_old * (V * si::cubic_metres / V_old) 
                    + dt * common::henry::mass_trans(rw2, D, acc_coeff, T, M_gas) 
                       * c * M_aq / M_gas * V * si::cubic_metres * rhod ) 
                   /
                   (real_t(1.) + common::henry::mass_trans(rw2, D, acc_coeff, T, M_gas) * dt 
                                 / common::henry::H_temp(T, H, dHR) / common::moist_air::kaBoNA<real_t>() / T)
                 ) / si::kilograms;
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_henry(
      const real_t &dt,
      const bool &chem_sys_cls
    )
    {   
      using namespace common::henry;      // H-prefixed
      using namespace common::molar_mass; // M-prefixed

      thrust_device::vector<real_t> &V(tmp_device_real_part);
      thrust_device::vector<real_t> &V_old(tmp_device_real_part_V_old);

      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      // gas absorption
      // TODO: so far it works only for 0D parcel model
      // TODO: open/close system logic -> async/sync timestep
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
      //molar mass of gases
      static const quantity<common::mass_over_amount, real_t> M_gas_[chem_gas_n] = {
        M_HNO3<real_t>(), M_NH3<real_t>(),  M_CO2<real_t>(),
        M_SO2<real_t>(),  M_H2O2<real_t>(), M_O3<real_t>()
      };
      //molar mass of dissolved chem species
      static const quantity<common::mass_over_amount, real_t> M_aq_[chem_gas_n] = {
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
          pi_t,                                              // trace gases
          typename thrust_device::vector<real_t>::iterator,  // m_old
          typename thrust_device::vector<real_t>::iterator,  // rw2
          pi_t,                                              // rhod
          typename thrust_device::vector<real_t>::iterator   // V old
        >
      > zip_it_t;

      if (chem_sys_cls == true){   //closed chemical system - reduce mixing ratio due to Henrys law

        // temporarily needed to store old mass per cell 
        thrust_device::vector<real_t> &tmp_mass(tmp_device_real_part_mass);
        std::map<enum chem_species_t, thrust_device::vector<real_t> > mass_old;
        std::map<enum chem_species_t, thrust_device::vector<real_t> > mass_new;

        for (int i = 0; i < chem_gas_n; ++i)
        {
          mass_old[(chem_species_t)i].resize(n_cell);
          mass_new[(chem_species_t)i].resize(n_cell);

          // store the total mass of chem species in cloud droplets per cell
          thrust::transform(
            n.begin(), n.end(),             // input - 1st arg
            chem_bgn[i],                    // input - 2nd arg
            tmp_mass.begin(),               // output
            detail::chem_summator<n_t, real_t>() // op
          );
          mass_old[(chem_species_t)i][0] = thrust::reduce(tmp_mass.begin(), tmp_mass.end(),(real_t) 0, thrust::plus<real_t>());

          // apply Henrys law tp the in-drop chemical compounds 
          thrust::transform(
            V.begin(), V.end(),            // input - 1st arg
            zip_it_t(thrust::make_tuple(   // input - 2nd arg
              thrust::make_permutation_iterator(p.begin(), ijk.begin()),
              thrust::make_permutation_iterator(T.begin(), ijk.begin()),
              thrust::make_permutation_iterator(ambient_chem[(chem_species_t)i].begin(), ijk.begin()),
              chem_bgn[i],
              rw2.begin(),
              thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
              V_old.begin()
            )),
            chem_bgn[i],                                                                                        // output
            detail::chem_Henry_fun<real_t>(H_[i], dHR_[i], M_gas_[i], M_aq_[i], D_[i], ac_[i], dt * si::seconds) // op
            );

          // store the total mass of chem species in cloud droplets per cell after Henry
          thrust::transform(
            n.begin(), n.end(),             // input - 1st arg
            chem_bgn[i],                    // input - 2nd arg
            tmp_mass.begin(),               // output
            detail::chem_summator<n_t, real_t>() // op
          );
          mass_new[(chem_species_t)i][0] = thrust::reduce(tmp_mass.begin(), tmp_mass.end(),(real_t) 0, thrust::plus<real_t>());

          // apply the change to the mixing ratios of trace gases
          ambient_chem[(chem_species_t)i][0] -= 
            (mass_new[(chem_species_t)i][0] - mass_old[(chem_species_t)i][0]) / M_aq_[i] * M_gas_[i] / dv[0] / rhod[0];
        }
      }
      else{ // open chemical system - do not change trace gase mixing ratios due to Henrys law
        for (int i = 0; i < chem_gas_n; ++i)
        {
          // apply Henrys law tp the in-drop chemical compounds 
          thrust::transform(
            V.begin(), V.end(),            // input - 1st arg
            zip_it_t(thrust::make_tuple(   // input - 2nd arg
              thrust::make_permutation_iterator(p.begin(), ijk.begin()),
              thrust::make_permutation_iterator(T.begin(), ijk.begin()),
              thrust::make_permutation_iterator(ambient_chem[(chem_species_t)i].begin(), ijk.begin()),
              chem_bgn[i],
              rw2.begin(),
              thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
              V_old.begin()
            )),
            chem_bgn[i],                                                                                        // output
            detail::chem_Henry_fun<real_t>(H_[i], dHR_[i], M_gas_[i], M_aq_[i], D_[i], ac_[i], dt * si::seconds) // op
          );
        }
      }
    }
  };  
};
