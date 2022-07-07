// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/henry.hpp>
#include <libcloudph++/common/dissoc.hpp>

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
      { // calculate the change in trace gases due to Henrys law
        const real_t M_gas; //quantity<common::mass_over_amount, real_t>
        const real_t M_aq;  //quantity<common::mass_over_amount, real_t> 
 
        // ctor
        ambient_chem_calculator(
          const real_t &M_aq,
          const real_t &M_gas
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

          quantity<si::dimensionless, real_t> new_c  = c - (m_new * si::kilograms - m_old) / M_aq * M_gas / dv / rhod;

          // As of now the aq. chemistry module computes Henrys law from the point of view of each super droplet.
          // It never checks if the sum of decrements of the ambient air trace gas concentarions
          // (the sum over all super droplets in each grid-box) does not exceede the total concentration in each drid-box.
          // It can only happen when ambient air concentrations are approaching zero, 
          // so the error introduced here is small. 
          // TODO - still it would be good to have some check build in Henrys law, instead of just an assert.
          return new_c > 0 ? new_c : new_c * real_t(0);
        }
      };

      template <typename real_t>
      struct chem_Henry_fun
      { // gas absorption into cloud droplets (Henrys law)
        const int chem_iter;
        const real_t H;         // quantity<common::amount_over_volume_over_pressure, real_t>
        const real_t dHR;       // quantity<si::temperature, real_t>
        const real_t M_gas;     // quantity<common::mass_over_amount, real_t>
        const real_t M_aq;      // quantity<common::mass_over_amount, real_t>
        const real_t D;         // quantity<common::diffusivity, real_t>
        const real_t acc_coeff; // quantity<si::dimensionless, real_t>
        const real_t dt;        // quantity<si::time, real_t>

        // ctor
        BOOST_GPU_ENABLED
        chem_Henry_fun(
          const int chem_iter,
          const real_t &H,
          const real_t &dHR,
          const real_t &M_gas,
          const real_t &M_aq,
          const real_t &D,
          const real_t &acc_coeff,
          const real_t &dt
        ) : 
          chem_iter(chem_iter), H(H), dHR(dHR), M_gas(M_gas), M_aq(M_aq), D(D), acc_coeff(acc_coeff), dt(dt)
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
          const quantity<si::mass, real_t>          m_H    = thrust::get<6>(tpl) * si::kilograms;     

          using namespace common::henry;      // H-prefixed
          using namespace common::molar_mass; // M-prefixed
          using namespace common::dissoc;     // K-prefixed

          //helper for mass of the un-dissociated chem. compounds
          quantity<si::dimensionless, real_t> hlp;

          // helper for H+ concentration
          quantity<common::amount_over_volume, real_t> conc_H;
          conc_H = m_H / M_H<real_t>() / (V * si::cubic_metres);

          // helper for Henry coefficient corrected for temperature and pH
          typedef divide_typeof_helper<
            divide_typeof_helper<
              si::amount,
              si::volume
            >::type,
            si::pressure
          >::type amount_over_volume_over_pressure;

          quantity<amount_over_volume_over_pressure, real_t> Henry;

          switch(chem_iter)
          {
            case SO2:
            {
              quantity<common::amount_over_volume, real_t> Kt_SO2, Kt_HSO3;

              Kt_SO2  = K_temp(T, K_SO2<real_t>(),  dKR_SO2<real_t>());
              Kt_HSO3 = K_temp(T, K_HSO3<real_t>(), dKR_HSO3<real_t>());
          
              hlp = (real_t(1) + Kt_SO2/conc_H + Kt_SO2 * Kt_HSO3 / conc_H / conc_H);

              Henry = common::henry::H_temp(T, H * si::moles / si::cubic_metres / si::pascals, dHR * si::kelvins) * hlp;
            }
            break;

            case CO2:
            {
              quantity<common::amount_over_volume, real_t> Kt_CO2, Kt_HCO3;

              Kt_CO2  = K_temp(T, K_CO2<real_t>(),  dKR_CO2<real_t>());
              Kt_HCO3 = K_temp(T, K_HCO3<real_t>(), dKR_HCO3<real_t>());

              hlp = (real_t(1) + Kt_CO2/conc_H + Kt_CO2 * Kt_HCO3 / conc_H / conc_H);

              Henry = common::henry::H_temp(T, H * si::moles / si::cubic_metres / si::pascals, dHR * si::kelvins) * hlp;
            }
            break;

            case HNO3:
            {
              quantity<common::amount_over_volume, real_t> Kt_HNO3;
              Kt_HNO3 = K_temp(T, K_HNO3<real_t>(), dKR_HNO3<real_t>());

              hlp = (real_t(1) + Kt_HNO3 / conc_H);

              Henry = common::henry::H_temp(T, H * si::moles / si::cubic_metres / si::pascals, dHR * si::kelvins) * hlp;
            }
            break;

            case NH3:
            {
              quantity<common::amount_over_volume, real_t> Kt_NH3;
              Kt_NH3  = K_temp(T, K_NH3<real_t>(),  dKR_NH3<real_t>());

              hlp = (real_t(1.) + Kt_NH3 / K_H2O<real_t>() * conc_H);

              Henry = common::henry::H_temp(T, H * si::moles / si::cubic_metres / si::pascals, dHR * si::kelvins) * hlp;
            }
            break;

            case O3:
            {
              Henry = common::henry::H_temp(T, H * si::moles / si::cubic_metres / si::pascals, dHR * si::kelvins);
            }
            break;

            case H2O2:
            {
              Henry = common::henry::H_temp(T, H * si::moles / si::cubic_metres / si::pascals, dHR * si::kelvins);
            }
            break;

            default:
              assert(false);
          }
        
          // implicit solution to the eq. 8.22 from chapter 8.4.2 
          // in Peter Warneck Chemistry of the Natural Atmosphere  
          real_t mass_helper =     
          (
            ( m_old 
                + 
                (dt * si::seconds) * (V * si::cubic_metres) 
                * common::henry::mass_trans(
                                              rw2, 
                                              (D * si::metres * si::metres / si::seconds), 
                                              (acc_coeff * si::seconds/si::seconds), 
                                              T, 
                                              (M_gas * si::kilograms / si::moles)
                                            ) 
                * c * rhod * (M_aq / M_gas)
              )
              /
              (real_t(1.) + (dt * si::seconds) 
                 * common::henry::mass_trans(rw2, (D * si::metres * si::metres / si::seconds), 
                                            acc_coeff * si::seconds / si::seconds, T, (M_gas * si::kilograms / si::moles)) 
                 / Henry / common::moist_air::kaBoNA<real_t>() / T)
          ) / si::kilograms ;

          return mass_helper;// > 0 ? mass_helper : 0;
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_henry(
      const real_t &dt
    )
    {   
      using namespace common::henry;      // H-prefixed
      using namespace common::molar_mass; // M-prefixed
      using namespace common::dissoc;     // K-prefixed

      const thrust_device::vector<unsigned int> &chem_flag(tmp_device_n_part);
      const thrust_device::vector<real_t> &V(tmp_device_real_part);

      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      // gas absorption
      assert(
        HNO3 == 0 && NH3  == 1 && CO2 == 2 &&
        SO2  == 3 && H2O2 == 4 && O3  == 5 
      );
      // Henry constant 
      // if Boost.units would work for Thrust: const quantity<common::amount_over_volume_over_pressure, real_t>
      const real_t H_[chem_gas_n] = {
        H_HNO3<real_t>() / si::moles * si::pascals * si::cubic_metres, 
        H_NH3<real_t>()  / si::moles * si::pascals * si::cubic_metres,  
        H_CO2<real_t>()  / si::moles * si::pascals * si::cubic_metres,
        H_SO2<real_t>()  / si::moles * si::pascals * si::cubic_metres,  
        H_H2O2<real_t>() / si::moles * si::pascals * si::cubic_metres, 
        H_O3<real_t>()   / si::moles * si::pascals * si::cubic_metres
      } ;
      // correction to Henry const due to temperature
      const real_t dHR_[chem_gas_n] = {
        dHR_HNO3<real_t>() / si::kelvins, 
        dHR_NH3<real_t>()  / si::kelvins,  
        dHR_CO2<real_t>()  / si::kelvins,
        dHR_SO2<real_t>()  / si::kelvins,  
        dHR_H2O2<real_t>() / si::kelvins, 
        dHR_O3<real_t>()   / si::kelvins
      };
      //molar mass of gases
      const real_t M_gas_[chem_gas_n] = {
        M_HNO3<real_t>() * si::moles / si::kilograms, 
        M_NH3<real_t>()  * si::moles / si::kilograms,  
        M_CO2<real_t>()  * si::moles / si::kilograms,
        M_SO2<real_t>()  * si::moles / si::kilograms,  
        M_H2O2<real_t>() * si::moles / si::kilograms, 
        M_O3<real_t>()   * si::moles / si::kilograms
      };
      //molar mass of dissolved chem species
      const real_t M_aq_[chem_gas_n] = {
        M_HNO3<real_t>()    * si::moles / si::kilograms,    
        M_NH3_H2O<real_t>() * si::moles / si::kilograms, 
        M_CO2_H2O<real_t>() * si::moles / si::kilograms,
        M_SO2_H2O<real_t>() * si::moles / si::kilograms, 
        M_H2O2<real_t>()    * si::moles / si::kilograms, 
        M_O3<real_t>()      * si::moles / si::kilograms
      };
      //gas phase diffusivity
      const real_t D_[chem_gas_n] = {
        D_HNO3<real_t>() * si::seconds / si::metres / si::metres, 
        D_NH3<real_t>()  * si::seconds / si::metres / si::metres,  
        D_CO2<real_t>()  * si::seconds / si::metres / si::metres,
        D_SO2<real_t>()  * si::seconds / si::metres / si::metres,  
        D_H2O2<real_t>() * si::seconds / si::metres / si::metres, 
        D_O3<real_t>()   * si::seconds / si::metres / si::metres
      };
      // accomodation coefficient
      const real_t ac_[chem_gas_n] = {
        ac_HNO3<real_t>(), 
        ac_NH3<real_t>(),  
        ac_CO2<real_t>(),
        ac_SO2<real_t>(),  
        ac_H2O2<real_t>(), 
        ac_O3<real_t>()
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

      //closed chemical system - reduce mixing ratio due to Henrys law
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
        thrust::transform_if(
          V.begin(), V.end(),                             // input - 1st arg
          thrust::make_zip_iterator(thrust::make_tuple(   // input - 2nd arg
            thrust::make_permutation_iterator(p.begin(), ijk.begin()),
            thrust::make_permutation_iterator(T_ref.begin(), ijk_ref_hlpr.begin()),
            thrust::make_permutation_iterator(ambient_chem[(chem_species_t)i].begin(), ijk.begin()),
            chem_bgn[i],
            rw2.begin(),
            thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
            chem_bgn[H]
          )),
          chem_flag.begin(),                                                                         //stencil
          chem_bgn[i],                                                                               // output
          detail::chem_Henry_fun<real_t>(i, H_[i], dHR_[i], M_gas_[i], M_aq_[i], D_[i], ac_[i], dt), // op
          thrust::identity<unsigned int>()
        );

        // store the total mass of chem species in cloud droplets per cell after Henry
        thrust::pair<
          typename thrust_device::vector<thrust_size_t>::iterator,
          typename thrust_device::vector<real_t>::iterator
        > it_pair =
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
        count_n = it_pair.first - count_ijk.begin();
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
    
#if !defined(__NVCC__)
      using boost::math::isfinite;
#endif
      for (int i = 0; i < chem_gas_n; ++i){
        //debug::print(chem_bgn[i], chem_end[i]);
        assert(isfinite(*thrust::min_element(chem_bgn[i], chem_end[i])));
      }
    }
  };  
};
