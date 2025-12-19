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
      struct is_positive //: public thrust::unary_function<real_t, real_t>
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
      struct precip_rate_ice
      {
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t> &tpl)  // tpl is a tuple (a, c, rho, vt)
        {
          return real_t(4./3)
          #if !defined(__NVCC__)
            * pi<real_t>()
          #else
            * CUDART_PI
          #endif
          * thrust::get<0>(tpl) * thrust::get<0>(tpl) * thrust::get<1>(tpl) // a^2 * c
          * thrust::get<2>(tpl)  // rho
          * thrust::get<3>(tpl); // vt
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
      struct get_sqrt// : public thrust::unary_function<real_t, real_t>
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

      template<class real_t>
      class ice_mass
      {
      public:
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl)  // tpl is a tuple (a, c, rho)
        {
          return real_t(4./3)
          #if !defined(__NVCC__)
            * pi<real_t>()
          #else
            * CUDART_PI
          #endif
          * thrust::get<0>(tpl) * thrust::get<0>(tpl) * thrust::get<1>(tpl) // a^2 * c
          * thrust::get<2>(tpl);  // rho
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

      thrust_device::vector<real_t> &n_filtered = pimpl->n_filtered_gp->get();

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
#if !defined(__NVCC__)
      using std::pow;
#endif
      pimpl->moms_rng(pow(r_min, 3), pow(r_max, 3), pimpl->rd3.begin(), false);
    }

    // selects particles with (r_w >= r_min && r_w < r_max)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_wet_rng(const real_t &r_min, const real_t &r_max)
    {
#if !defined(__NVCC__)
      using std::pow;
#endif
      pimpl->moms_rng(pow(r_min, 2), pow(r_max, 2), pimpl->rw2.begin(), false);
    }

    // selects particles with (ice_a >= a_min && ice_a < a_max)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_ice_a_rng(const real_t &a_min, const real_t &a_max)
    {
      pimpl->moms_rng(a_min, a_max, pimpl->ice_a.begin(), false);
    }

    // selects particles with (ice_c >= c_min && ice_c < c_max)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_ice_c_rng(const real_t &c_min, const real_t &c_max)
    {
      pimpl->moms_rng(c_min, c_max, pimpl->ice_c.begin(), false);
    }

    // selects particles with (kpa >= kpa_min && kpa < kpa_max)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_kappa_rng(const real_t &kpa_min, const real_t &kpa_max)
    {
      pimpl->moms_rng(kpa_min, kpa_max, pimpl->kpa.begin(), false);
    }

    // selects ice particles 
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_ice()
    {
      pimpl->moms_gt0(pimpl->ice_a.begin()); // ice_a greater than 0
    }

    // selects water particles 
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_water()
    {
      pimpl->moms_eq0(pimpl->ice_a.begin()); // ice_a equal to 0
    }

    // selects particles with (r_d >= r_min && r_d < r_max) from particles previously selected
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_dry_rng_cons(const real_t &r_min, const real_t &r_max)
    {
#if !defined(__NVCC__)
      using std::pow;
#endif
      pimpl->moms_rng(pow(r_min, 3), pow(r_max, 3), pimpl->rd3.begin(), true);
    }

    // selects particles with (r_w >= r_min && r_w < r_max) from particles previously selected
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_wet_rng_cons(const real_t &r_min, const real_t &r_max)
    {
#if !defined(__NVCC__)
      using std::pow;
#endif
      pimpl->moms_rng(pow(r_min, 2), pow(r_max, 2), pimpl->rw2.begin(), true);
    }

    // selects particles with (ice_a >= a_min && ice_a < a_max) from particles previously selected
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_ice_a_rng_cons(const real_t &a_min, const real_t &a_max)
    {
      pimpl->moms_rng(a_min, a_max, pimpl->ice_a.begin(), true);
    }

    // selects particles with (ice_c >= c_min && ice_c < c_max) from particles previously selected
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_ice_c_rng_cons(const real_t &c_min, const real_t &c_max)
    {
      pimpl->moms_rng(c_min, c_max, pimpl->ice_c.begin(), true);
    }

    // selects particles with (kpa >= kpa_min && kpa < kpa_max) from particles previously selected
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_kappa_rng_cons(const real_t &kpa_min, const real_t &kpa_max)
    {
      pimpl->moms_rng(kpa_min, kpa_max, pimpl->kpa.begin(), true);
    }

    // selects ice particles from particles previously selected
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_ice_cons()
    {
      pimpl->moms_gt0(pimpl->ice_a.begin(), true); // ice_a greater than 0
    }

    // selects water particles from particles previously selected 
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_water_cons()
    {
      pimpl->moms_eq0(pimpl->ice_a.begin(), true); // ice_a equal to 0
    }

    // selects particles with RH >= Sc   (Sc - critical supersaturation)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_RH_ge_Sc()
    {
      // intentionally using the same tmp vector as inside moms_cmp below
      // thrust_device::vector<real_t> &RH_minus_Sc(pimpl->tmp_device_real_part);
      auto RH_minus_Sc_g = pimpl->tmp_device_real_part.get_guard();
      thrust_device::vector<real_t> &RH_minus_Sc = RH_minus_Sc_g.get();

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
      // thrust_device::vector<real_t> &rc2(pimpl->tmp_device_real_part);
      auto rc2_g = pimpl->tmp_device_real_part.get_guard();
      thrust_device::vector<real_t> &rc2 = rc2_g.get();

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

    // computes n-th moment of the ice equatorial radius spectrum for the selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_ice_a_mom(const int &n)
    {
      pimpl->moms_calc(pimpl->ice_a.begin(), n);
    }

    // computes n-th moment of the ice polar radius spectrum for the selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_ice_c_mom(const int &n)
    {
      pimpl->moms_calc(pimpl->ice_c.begin(), n);
    }

     // computes ice mixing ratio
    template <typename real_t, backend_t device>
      void particles_t<real_t, device>::diag_ice_mix_ratio()
     {
      pimpl->moms_calc(thrust::make_transform_iterator(
        thrust::make_zip_iterator(thrust::make_tuple(pimpl->ice_a.begin(), pimpl->ice_c.begin(), pimpl->ice_rho.begin())),
        detail::ice_mass<real_t>()
        ),
        real_t(1));
     }

    // compute n-th moment of kappa for selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_kappa_mom(const int &n)
    {   
      pimpl->moms_calc(pimpl->kpa.begin(), n);
    }   

    // compute n-th moment of SGS velocity in x for selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_up_mom(const int &n)
    {   
      pimpl->moms_calc(pimpl->up.begin(), n);
    }   

    // compute n-th moment of SGS velocity in y for selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_vp_mom(const int &n)
    {   
      pimpl->moms_calc(pimpl->vp.begin(), n);
    }   

    // compute n-th moment of SGS velocity in z for selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_wp_mom(const int &n)
    {   
      pimpl->moms_calc(pimpl->wp.begin(), n);
    }   

    // compute n-th moment of incloud_time for selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_incloud_time_mom(const int &n)
    {   
      if(pimpl->opts_init.diag_incloud_time)
        pimpl->moms_calc(pimpl->incloud_time.begin(), n);
      else
        assert(0 && "diag_incloud_time_mom called, but opts_init.diag_incloud_time==false");
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
      auto tmp_vt_g = pimpl->tmp_device_real_part.get_guard();
      thrust_device::vector<real_t> &tmp_vt = tmp_vt_g.get();
      thrust::copy(pimpl->vt.begin(), pimpl->vt.end(), tmp_vt.begin());
    
      thrust::transform(
        pimpl->vt.begin(),
        pimpl->vt.end(),
        pimpl->rw2.begin(),
        pimpl->vt.begin(),
        detail::precip_rate<real_t>()
      );  

      pimpl->moms_all(); // we need this here, because hskpng_vterm modifies tmp_device_real_part, which is used as n_filtered in moms_calc
      pimpl->moms_calc(pimpl->vt.begin(), 1., false);
 
      // copy back stored vterm
      thrust::copy(tmp_vt.begin(), tmp_vt.end(), pimpl->vt.begin());
    }

    // compute 1st (non-specific) moment of ice_mass * vt of all SDs
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_precip_rate_ice_mass()
    {
      // updating terminal velocities
      pimpl->hskpng_vterm_all();

      pimpl->moms_all(); // we need this here, because hskpng_vterm modifies tmp_device_real_part, which is used as n_filtered in moms_calc
      pimpl->moms_calc(thrust::make_transform_iterator(
            thrust::make_zip_iterator(thrust::make_tuple(pimpl->ice_a.begin(), pimpl->ice_c.begin(),
              pimpl->ice_rho.begin(), pimpl->vt.begin())),
            detail::precip_rate_ice<real_t>()
            ),
            real_t(1), false);
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
      if(pimpl->opts_init.chem_switch == false) throw std::runtime_error("libcloudph++: all chemistry was switched off in opts_init");
      pimpl->moms_calc(pimpl->chem_bgn[c], 1.);
    }

    template <typename real_t, backend_t device>
    std::map<common::output_t, real_t> particles_t<real_t, device>::diag_puddle()
    {
      return pimpl->output_puddle;
    }
  };
};
