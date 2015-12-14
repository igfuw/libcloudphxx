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
      template <typename real_t>
      struct chem_summator
      { // calculate the mass of chem compounds (multiplicity * mass)  
        template <typename tup_t>
        BOOST_GPU_ENABLED
        real_t operator()(const tup_t &tpl) const
        {
          return thrust::get<0>(tpl) * thrust::get<1>(tpl);
        }
      };

      template <typename real_t>
      struct ambient_chem_calculator
      { // calculate the chenge in trace gases due to Henrys law
        const quantity<common::mass_over_amount, real_t> M_gas;
        const quantity<common::mass_over_amount, real_t> M_aq;
 
        // ctor
        ambient_chem_calculator(
          const quantity<common::mass_over_amount, real_t> &M_aq,
          const quantity<common::mass_over_amount, real_t> &M_gas
        ) : 
          M_gas(M_gas), M_aq(M_aq)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(
          const real_t &m_new, 
          const thrust::tuple<real_t, real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::mass, real_t>          m_old  = thrust::get<0>(tpl) * si::kilograms;     
          const quantity<si::mass_density, real_t>  rhod   = thrust::get<1>(tpl) * si::kilograms / si::cubic_metres; 
          const quantity<si::volume, real_t>        dv     = thrust::get<2>(tpl) * si::cubic_metres;
          const quantity<si::dimensionless, real_t> c      = thrust::get<3>(tpl); 

          return c - (m_new * si::kilograms - m_old) / M_aq * M_gas / dv / rhod;
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

          // implicit solution to the eq. 8.22 from chapter 8.4.2 
          // in Peter Warneck Chemistry of the Natural Atmosphere  
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

      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      // gas absorption
      // TODO: open/close system logic -> async/sync timestep
      assert(
        HNO3 == 0 && NH3  == 1 && CO2 == 2 &&
        SO2  == 3 && H2O2 == 4 && O3  == 5 
      );
      //Henry constant
      const quantity<common::amount_over_volume_over_pressure, real_t> H_[chem_gas_n] = {
        H_HNO3<real_t>(), H_NH3<real_t>(),  H_CO2<real_t>(),
        H_SO2<real_t>(),  H_H2O2<real_t>(), H_O3<real_t>()
      };
      //correction to Henry const due to temperature
      const quantity<si::temperature, real_t> dHR_[chem_gas_n] = {
        dHR_HNO3<real_t>(), dHR_NH3<real_t>(),  dHR_CO2<real_t>(),
        dHR_SO2<real_t>(),  dHR_H2O2<real_t>(), dHR_O3<real_t>()
      };
      //molar mass of gases
      const quantity<common::mass_over_amount, real_t> M_gas_[chem_gas_n] = {
        M_HNO3<real_t>(), M_NH3<real_t>(),  M_CO2<real_t>(),
        M_SO2<real_t>(),  M_H2O2<real_t>(), M_O3<real_t>()
      };
      //molar mass of dissolved chem species
      const quantity<common::mass_over_amount, real_t> M_aq_[chem_gas_n] = {
        M_HNO3<real_t>(),    M_NH3_H2O<real_t>(), M_CO2_H2O<real_t>(),
        M_SO2_H2O<real_t>(), M_H2O2<real_t>(), M_O3<real_t>()
      };
      //gas phase diffusivity
      const quantity<common::diffusivity, real_t> D_[chem_gas_n] = {
        D_HNO3<real_t>(), D_NH3<real_t>(),  D_CO2<real_t>(),
        D_SO2<real_t>(),  D_H2O2<real_t>(), D_O3<real_t>()
      };
      //accomodation coefficient
      const quantity<si::dimensionless, real_t> ac_[chem_gas_n] = {
        ac_HNO3<real_t>(), ac_NH3<real_t>(),  ac_CO2<real_t>(),
        ac_SO2<real_t>(),  ac_H2O2<real_t>(), ac_O3<real_t>()
      };

      // helpers for sorting out droplets
      typedef thrust::permutation_iterator<
        typename thrust_device::vector<n_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_n_t;
      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_chem_t;
      typedef thrust::zip_iterator<thrust::tuple<pi_n_t, pi_chem_t> > zip_it_t;

      if (chem_sys_cls == true){   //closed chemical system - reduce mixing ratio due to Henrys law

        hskpng_sort();

        // temporarily needed to store old mass per cell 
        thrust_device::vector<real_t> &mass_old(tmp_device_real_cell);
        thrust_device::vector<real_t> &mass_new(tmp_device_real_cell1);

        for (int i = 0; i < chem_gas_n; ++i)
        {
          // store the total mass of chem species in cloud droplets per cell
          thrust::reduce_by_key(
            sorted_ijk.begin(), sorted_ijk.end(),
            thrust::transform_iterator<             
              detail::chem_summator<real_t>,
              zip_it_t,
              real_t
            >(
              zip_it_t(thrust::make_tuple(
                pi_n_t(n.begin(), sorted_id.begin()),
                pi_chem_t(chem_bgn[i], sorted_id.begin())
              )),
              detail::chem_summator<real_t>() 
            ),
            count_ijk.begin(),
            mass_old.begin()
          );

          // apply Henrys law to the in-drop chemical compounds 
          thrust::transform(
            V.begin(), V.end(),                             // input - 1st arg
            thrust::make_zip_iterator(thrust::make_tuple(   // input - 2nd arg
              thrust::make_permutation_iterator(p.begin(), ijk.begin()),
              thrust::make_permutation_iterator(T.begin(), ijk.begin()),
              thrust::make_permutation_iterator(ambient_chem[(chem_species_t)i].begin(), ijk.begin()),
              chem_bgn[i],
              rw2.begin(),
              thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
              V_old.begin()
            )),
            chem_bgn[i],                                                                                         // output
            detail::chem_Henry_fun<real_t>(H_[i], dHR_[i], M_gas_[i], M_aq_[i], D_[i], ac_[i], dt * si::seconds) // op
          );

          //TODO - check if its needed in some 2D test
          /*thrust::transform(
            thrust::make_permutation_iterator(V.begin(), sorted_id.begin()),
            thrust::make_permutation_iterator(V.begin(), sorted_id.end()),
            thrust::make_zip_iterator(thrust::make_tuple(   // input - 2nd arg
              thrust::make_permutation_iterator(p.begin(), sorted_ijk.begin()),
              thrust::make_permutation_iterator(T.begin(), sorted_ijk.begin()),
              thrust::make_permutation_iterator(ambient_chem[(chem_species_t)i].begin(), sorted_ijk.begin()),
              thrust::make_permutation_iterator(chem_bgn[i].begin(), sorted_id.begin()),
              thrust::make_permutation_iterator(rw2.begin(), sorted_id.begin()),
              thrust::make_permutation_iterator(rhod.begin(), sorted_ijk.begin()),
              thrust::make_permutation_iterator(V_old.begin(), sorted_id.begin())
            )),
            thrust::make_permutation_iterator(chem_bgn[i], sorted_id.begin()),
            detail::chem_Henry_fun<real_t>(H_[i], dHR_[i], M_gas_[i], M_aq_[i], D_[i], ac_[i], dt * si::seconds) // op
          );*/

          // store the total mass of chem species in cloud droplets per cell after Henry
          thrust::pair<
            typename thrust_device::vector<thrust_size_t>::iterator,
            typename thrust_device::vector<real_t>::iterator
          > np =
          thrust::reduce_by_key(
            sorted_ijk.begin(), sorted_ijk.end(),
            thrust::transform_iterator<             // input - values
              detail::chem_summator<real_t>,
              zip_it_t,
              real_t
            >(
              zip_it_t(thrust::make_tuple(
                pi_n_t(n.begin(), sorted_id.begin()),
                pi_chem_t(chem_bgn[i], sorted_id.begin())
              )),
              detail::chem_summator<real_t>()      // op
            ),
            count_ijk.begin(),
            mass_new.begin()
          );
          count_n = np.first - count_ijk.begin();
          assert(count_n > 0 && count_n <= n_cell);


          // apply the change to the mixing ratios of trace gases
          thrust::transform(
            mass_new.begin(), mass_new.begin() + count_n,                              // input - 1st arg
            thrust::make_zip_iterator(thrust::make_tuple(                              // input - 2nd arg
              mass_old.begin(),
              thrust::make_permutation_iterator(rhod.begin(), count_ijk.begin()), 
              thrust::make_permutation_iterator(dv.begin(), count_ijk.begin()),
              thrust::make_permutation_iterator(ambient_chem[(chem_species_t)i].begin(), count_ijk.begin())
            )),
            thrust::make_permutation_iterator(ambient_chem[(chem_species_t)i].begin(), count_ijk.begin()), // output 
            detail::ambient_chem_calculator<real_t>(M_aq_[i], M_gas_[i])                                   // op
          );
          assert(*thrust::min_element(
            ambient_chem[(chem_species_t)i].begin(), ambient_chem[(chem_species_t)i].end()
          ) >= 0);
        }
      }

      else{ // open chemical system - do not change trace gase mixing ratios due to Henrys law
        for (int i = 0; i < chem_gas_n; ++i)
        {
          // apply Henrys law to the in-drop chemical compounds 
          thrust::transform(
            V.begin(), V.end(),                             // input - 1st arg
            thrust::make_zip_iterator(thrust::make_tuple(   // input - 2nd arg
              thrust::make_permutation_iterator(p.begin(), ijk.begin()),
              thrust::make_permutation_iterator(T.begin(), ijk.begin()),
              thrust::make_permutation_iterator(ambient_chem[(chem_species_t)i].begin(), ijk.begin()),
              chem_bgn[i],
              rw2.begin(),
              thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
              V_old.begin()
            )),
            chem_bgn[i],                                                                                         // output
            detail::chem_Henry_fun<real_t>(H_[i], dHR_[i], M_gas_[i], M_aq_[i], D_[i], ac_[i], dt * si::seconds) // op
          );
        }
      }

    }
  };  
};
