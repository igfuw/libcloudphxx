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
      struct adj_rw2
      {
        template<class real_t>
        BOOST_GPU_ENABLED
        bool operator()(real_t rw2, real_t m_added)
        {
          // turn into added volume
          m_added /= libcloudphxx::common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres;
          return pow(pow(rw2, real_t(3./2)) + 3. / 4. / 
#if !defined(__NVCC__)
            pi<real_t>()
#else
            CUDART_PI
#endif
            * m_added, real_t(2./3));
        }
      };

      struct is_greater_than_zero
      {
        template<class real_t>
        BOOST_GPU_ENABLED
        bool operator()(real_t x)
        {
          return x > real_t(0.);
        }
      };
      struct divide
      {
        template<class real_t, class n_t>
        BOOST_GPU_ENABLED
        bool operator()(real_t x, n_t y)
        {
          return x / real_t(y);
        }
      };
    };
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::rc_adjust(
      const arrinfo_t<real_t> &_rc_adjust
    ) 
    {   
      // use tmp vector to store adjustment
      thrust_device::vector<real_t> &rc_adjust(tmp_device_real_cell);

      // save input adjustment to the temp vector
      sync(_rc_adjust, rc_adjust);

      // alias for filtered SDs with RH>=Sc
      thrust_device::vector<real_t> &n_filtered(tmp_device_real_part);

      // turn n_filtered into activated? = yes/no 
      thrust::replace_if(
        n_filtered.begin(),
        n_filtered.end(),
        detail::is_greater_than_zero(),
        real_t(1)
      );

      hskpng_sort();

      // count number of activated SDs per cell
      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        thrust_device::vector<n_t>::iterator
      > n = thrust::reduce_by_key(
        sorted_ijk.begin(), sorted_ijk.end(),   // input - keys
        n_filtered.begin(),                     // input - values
        count_ijk.begin(),                      // output - keys
        count_num.begin()                       // output - values
      );  
      count_n = n.first - count_ijk.begin();

      // divide water adjustment by number of activated SDs
      thrust::transform(
        thrust::make_permutation_iterator(rc_adjust.begin(), count_ijk.begin()),
        thrust::make_permutation_iterator(rc_adjust.begin(), count_ijk.end()),
        thrust::make_permutation_iterator(count_num.begin(), count_ijk.begin()),
        thrust::make_permutation_iterator(rc_adjust.begin(), count_ijk.begin()),
        detail::divide()
      );

      // adjustment is in mass_of_added_water[kg] / mass_of_dry_air[kg]
      // multiply it now by mass of dry air in the cell
      thrust::transform(
        rc_adjust.begin(), rc_adjust.end(),
        dv.begin(),
        rc_adjust.begin(),
        thrust::multiplies<real_t>()
      );
      thrust::transform(
        rc_adjust.begin(), rc_adjust.end(),
        rhod.begin(),
        rc_adjust.begin(),
        thrust::multiplies<real_t>()
      );

      // apply the adjustment - change rw of activated SDs
      thrust::transform(
        thrust::make_permutation_iterator(rw2.begin(), sorted_ijk.begin()),
        thrust::make_permutation_iterator(rw2.begin(), sorted_ijk.end()),
        thrust::make_permutation_iterator(rc_adjust.begin(), sorted_ijk.begin()), // mass to be added to each SD
        thrust::make_permutation_iterator(rw2.begin(), sorted_ijk.begin()),
        detail::adj_rw2()
      );
    }
  };  
};
