// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/kappa_koehler.hpp> // TODO: not here...

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail 
    {
      template <typename real_t>
      struct rw3_cr
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rd3, const thrust::tuple<real_t, real_t> &tpl)
        {
          const quantity<si::dimensionless, real_t> kpa = thrust::get<0>(tpl);
          const quantity<si::temperature, real_t> T = thrust::get<1>(tpl) * si::kelvins;

#if !defined(__NVCC__)
          using std::pow;
#endif
          return pow(
            common::kappa_koehler::rw3_cr(
              rd3 * si::cubic_metres, 
              kpa,
              T
            ) / si::cubic_metres
            , 
            real_t(2./3)
          ); 
        }
      };

      template <typename real_t>
      struct precip_rate
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &vt, const real_t &rw2)
        {
#if !defined(__NVCC__)
          using std::pow;
#endif
          return pow(rw2, real_t(3./2)) * vt;
        }
      };

      template <typename real_t>
      struct RH_minus_Sc
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rd3, const thrust::tuple<real_t, real_t, real_t> &tpl)
        {
          const quantity<si::dimensionless, real_t> kpa = thrust::get<0>(tpl);
          const quantity<si::temperature, real_t> T = thrust::get<1>(tpl) * si::kelvins;
          const quantity<si::dimensionless, real_t> RH = thrust::get<2>(tpl);

          return RH - common::kappa_koehler::S_cr(
	    rd3 * si::cubic_metres, 
	    kpa,
	    T
          );
        }
      };
    }

    // records relative humidity
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_RH()
    {
      pimpl->hskpng_Tpr(); 

      thrust::copy(
        pimpl->RH.begin(), 
        pimpl->RH.end(), 
        pimpl->count_mom.begin()
      );
    }

    // records super-droplet concentration per grid cell
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_sd_conc()
    {
      // common code with coalescence, hence separated into a method
      pimpl->hskpng_count(); 

      // n_t -> real_t cast
      thrust::copy(
        pimpl->count_num.begin(), 
        pimpl->count_num.end(), 
        pimpl->count_mom.begin()
      );
    }

    // selected all particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_all()
    {
      pimpl->moms_all();
    }

    // selects particles with (r_d >= r_min && r_d < r_max)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_dry_rng(const real_t &r_min, const real_t &r_max)
    {
      pimpl->moms_rng(pow(r_min, 3), pow(r_max, 3), pimpl->rd3.begin());
    }

    // selects particles with (r_w >= r_min && r_w < r_max)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_wet_rng(const real_t &r_min, const real_t &r_max)
    {
      pimpl->moms_rng(pow(r_min, 2), pow(r_max, 2), pimpl->rw2.begin());
    }

    // selects particles with RH >= Sc   (Sc - critical supersaturation)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_RH_ge_Sc()
    {
      // intentionally using the same tmp vector as inside moms_cmp below
      thrust_device::vector<real_t> &RH_minus_Sc(pimpl->tmp_device_real_part);

      // computing RH_minus_Sc for each particle
      thrust::transform(
        pimpl->rd3.begin(), pimpl->rd3.end(), // input - 1st arg
        thrust::make_zip_iterator(make_tuple(
          pimpl->kpa.begin(), 
          thrust::make_permutation_iterator(
            pimpl->T.begin(),
            pimpl->ijk.begin()
          ),
          thrust::make_permutation_iterator(
            pimpl->RH.begin(),
            pimpl->ijk.begin()
          )
        )),                                   // input - 2nd arg 
        RH_minus_Sc.begin(),                  // output
        detail::RH_minus_Sc<real_t>()         // op
      );

      // selecting those with RH - Sc >= 0
      pimpl->moms_ge0(RH_minus_Sc.begin());
    }

    // selects particles with rw >= rc   (rc - critical radius)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_rw_ge_rc()
    {
      // intentionally using the same tmp vector as inside moms_cmp below
      thrust_device::vector<real_t> &rc2(pimpl->tmp_device_real_part);

      // computing rc2 for each particle
      thrust::transform(
        pimpl->rd3.begin(), pimpl->rd3.end(), // input - 1st arg
        thrust::make_zip_iterator(make_tuple(
          pimpl->kpa.begin(), 
          thrust::make_permutation_iterator(
            pimpl->T.begin(),
            pimpl->ijk.begin()
          )
        )),                                   // input - 2nd arg 
        rc2.begin(),                          // output
        detail::rw3_cr<real_t>()              // op
      );

      // selecting those with rw2 >= rc2
      pimpl->moms_cmp(pimpl->rw2.begin(), rc2.begin());
    }

    // computes n-th moment of the dry spectrum for the selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_dry_mom(const int &n)
    {
      pimpl->moms_calc(pimpl->rd3.begin(), n/3.);
    }

    // computes n-th moment of the wet spectrum for the selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_wet_mom(const int &n)
    {
      pimpl->moms_calc(pimpl->rw2.begin(), n/2.);
    }

    // computes mass density function for wet radii using estimator from Shima et al. (2009)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_wet_mass_dens(const real_t &rad, const real_t &sig0)
    {
      pimpl->mass_dens_estim(pimpl->rw2.begin(), rad, sig0, 1./2.);
    }

    // compute 1st (non-specific) moment of rw^3 * vt
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_precip_rate()
    {   
      // updating terminal velocities
      pimpl->hskpng_vterm_all();

      // temporary vector to store vt
      thrust::host_vector<real_t> tmp_vt(pimpl->n_part);
      thrust::copy(pimpl->vt.begin(), pimpl->vt.end(), tmp_vt.begin());
    
      thrust::transform(
        pimpl->vt.begin(),
        pimpl->vt.end(),
        pimpl->rw2.begin(),
        pimpl->vt.begin(),
        detail::precip_rate<real_t>()
      );  

      pimpl->moms_calc_cond(pimpl->vt.begin(), 1.);
 
      // copy back stored vterm
      thrust::copy(tmp_vt.begin(), tmp_vt.end(), pimpl->vt.begin());
      // release the memory
      tmp_vt.erase(tmp_vt.begin(), tmp_vt.end());
    }   

    // computes mean chemical properties for the selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_chem(const enum chem_species_t &c)
    {
      if(pimpl->opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off in opts_init");
      pimpl->moms_calc(pimpl->chem_bgn[c], 1);
    }
  };
};
