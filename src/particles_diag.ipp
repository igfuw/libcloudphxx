// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/kappa_koehler.hpp> // TODO: not here...
#include <thrust/sequence.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail 
    {
      template <typename real_t>
      struct add_div
      {
        const real_t dx;
        add_div(const real_t &dx) : dx(dx) {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &div, const thrust::tuple<real_t, real_t> &tpl)
        {
          return div + (thrust::get<1>(tpl) - thrust::get<0>(tpl)) / dx;
        }
      };

      template <typename real_t>
      struct is_positive : public thrust::unary_function<real_t, real_t>
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &a)
        {
          return a > 0. ? 1 : 0;
        }
      };

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

      template <typename real_t>
      struct get_sqrt : public thrust::unary_function<real_t, real_t>
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2)
        {
#if !defined(__NVCC__)
          using std::sqrt;
#endif
          return sqrt(rw2);
        }
      };
    }

    // records pressure
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_pressure()
    {
      pimpl->hskpng_Tpr(); 

      thrust::copy(
        pimpl->p.begin(), 
        pimpl->p.end(), 
        pimpl->count_mom.begin()
      );

      // p defined in all cells
      pimpl->count_n = pimpl->n_cell;
      thrust::sequence(pimpl->count_ijk.begin(), pimpl->count_ijk.end());
    }

    // records temperature
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_temperature()
    {
      pimpl->hskpng_Tpr(); 

      thrust::copy(
        pimpl->T.begin(), 
        pimpl->T.end(), 
        pimpl->count_mom.begin()
      );

      // T defined in all cells
      pimpl->count_n = pimpl->n_cell;
      thrust::sequence(pimpl->count_ijk.begin(), pimpl->count_ijk.end());
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

      // RH defined in all cells
      pimpl->count_n = pimpl->n_cell;
      thrust::sequence(pimpl->count_ijk.begin(), pimpl->count_ijk.end());
    }

    // records super-droplet concentration per grid cell
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_sd_conc()
    {
      namespace arg = thrust::placeholders;
      assert(pimpl->selected_before_counting);

      thrust_device::vector<real_t> &n_filtered(pimpl->tmp_device_real_part);

      // similar to hskpng_count
      pimpl->hskpng_sort();

      // computing count_* - number of particles per grid cell
      auto n = thrust::reduce_by_key(
        pimpl->sorted_ijk.begin(), pimpl->sorted_ijk.end(),   // input - keys
        thrust::make_permutation_iterator(
          thrust::make_transform_iterator(n_filtered.begin(), detail::is_positive<real_t>()),
          pimpl->sorted_id.begin()
        ),
        pimpl->count_ijk.begin(),                      // output - keys
        pimpl->count_mom.begin()                      // output - values
      );
      pimpl->count_n = n.first - pimpl->count_ijk.begin();
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

    // selects particles with (kpa >= kpa_min && kpa < kpa_max)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_kappa_rng(const real_t &kpa_min, const real_t &kpa_max)
    {
      pimpl->moms_rng(kpa_min, kpa_max, pimpl->kpa.begin());
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

    // compute n-th moment of kappa for selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_kappa_mom(const int &n)
    {   
      pimpl->moms_calc(pimpl->kpa.begin(), n);
    }   

    // computes mass density function for wet radii using estimator from Shima et al. (2009)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_wet_mass_dens(const real_t &rad, const real_t &sig0)
    {
      pimpl->mass_dens_estim(pimpl->rw2.begin(), rad, sig0, 1./2.);
    }

    // to diagnose if velocity field is nondivergent
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_vel_div()
    {   
      if(pimpl->n_dims==0) return;

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<thrust_size_t>::iterator,
        typename thrust::counting_iterator<thrust_size_t>
      > pi;


      thrust::fill(
        pimpl->count_mom.begin(),
        pimpl->count_mom.end(),
        real_t(0.)
      );

      switch (pimpl->n_dims)
      {
        case 3:
          thrust::transform(
            pimpl->count_mom.begin(), //arg1
            pimpl->count_mom.begin() + pimpl->n_cell,
            thrust::make_zip_iterator(thrust::make_tuple(
              thrust::make_permutation_iterator(pimpl->courant_y.begin(), pi(pimpl->fre.begin(), thrust::make_counting_iterator<thrust_size_t>(pimpl->halo_x))),   // fre counts from the start of halo, but here we need only real cells
              thrust::make_permutation_iterator(pimpl->courant_y.begin(), pi(pimpl->hnd.begin(), thrust::make_counting_iterator<thrust_size_t>(pimpl->halo_x)))
            )), // arg2
            pimpl->count_mom.begin(), //out
            detail::add_div<real_t>(pimpl->opts_init.dt)
          );
        case 2:
          thrust::transform(
            pimpl->count_mom.begin(), //arg1
            pimpl->count_mom.begin() + pimpl->n_cell,
            thrust::make_zip_iterator(thrust::make_tuple(
              thrust::make_permutation_iterator(pimpl->courant_z.begin(), pi(pimpl->blw.begin(), thrust::make_counting_iterator<thrust_size_t>(pimpl->halo_x))),
              thrust::make_permutation_iterator(pimpl->courant_z.begin(), pi(pimpl->abv.begin(), thrust::make_counting_iterator<thrust_size_t>(pimpl->halo_x)))
            )), // arg2
            pimpl->count_mom.begin(), //out
            detail::add_div<real_t>(pimpl->opts_init.dt)
          );
        case 1:
          thrust::transform(
            pimpl->count_mom.begin(), //arg1
            pimpl->count_mom.begin() + pimpl->n_cell,
            thrust::make_zip_iterator(thrust::make_tuple(
              thrust::make_permutation_iterator(pimpl->courant_x.begin(), pi(pimpl->lft.begin(), thrust::make_counting_iterator<thrust_size_t>(pimpl->halo_x))),
              thrust::make_permutation_iterator(pimpl->courant_x.begin(), pi(pimpl->rgt.begin(), thrust::make_counting_iterator<thrust_size_t>(pimpl->halo_x)))
            )), // arg2
            pimpl->count_mom.begin(), //out
            detail::add_div<real_t>(pimpl->opts_init.dt)
          );
      }
      // divergence defined in all cells
      pimpl->count_n = pimpl->n_cell;
      thrust::sequence(pimpl->count_ijk.begin(), pimpl->count_ijk.end());
    }

    // compute 1st (non-specific) moment of rw^3 * vt of all SDs
    // TODO: replace it with simple diag vt?
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

      pimpl->moms_all();
      pimpl->moms_calc(pimpl->vt.begin(), 1., false);
 
      // copy back stored vterm
      thrust::copy(tmp_vt.begin(), tmp_vt.end(), pimpl->vt.begin());
      // release the memory
      tmp_vt.erase(tmp_vt.begin(), tmp_vt.end());
    }   

    // get max rw in each cell
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_max_rw()
    {   
      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::const_iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_t;

      pimpl->hskpng_sort();

      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        typename thrust_device::vector<real_t>::iterator
      > n = thrust::reduce_by_key(
        // input - keys
        pimpl->sorted_ijk.begin(), pimpl->sorted_ijk.end(),  
        // input - values
        thrust::make_transform_iterator(
          pi_t(pimpl->rw2.begin(),   pimpl->sorted_id.begin()),
          detail::get_sqrt<real_t>()
        ),
        // output - keys
        pimpl->count_ijk.begin(),
        // output - values
        pimpl->count_mom.begin(),
        // key comparison
        thrust::equal_to<real_t>(),
        // reduction type
        thrust::maximum<real_t>()
      );  

      pimpl->count_n = n.first - pimpl->count_ijk.begin();
      assert(pimpl->count_n > 0 && pimpl->count_n <= pimpl->n_cell);
    }

    // computes mean chemical properties for the selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_chem(const enum chem_species_t &c)
    {
      if(pimpl->opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off in opts_init");
      pimpl->moms_calc(pimpl->chem_bgn[c], 1.);
    }

    template <typename real_t, backend_t device>
    std::map<libcloudphxx::common::output_t, real_t> particles_t<real_t, device>::diag_puddle()
    {
      return pimpl->output_puddle;
    }
  };
};
