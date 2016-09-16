// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/reduce.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/zip_iterator.h>

namespace libcloudphxx
{
  namespace lgrngn
  {

    namespace detail
    {
      // selects particles with x >= min and x < max
      template <typename real_t>
      struct range_filter
      {
        real_t min, max;

        range_filter(real_t min, real_t max) : min(min), max(max) {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &y, const real_t &x)
        {
          return x >= min && x < max ? y : 0; 
        }
      };

      // selects particles with x >= c
      template <typename real_t>
      struct comp_filter
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &y, const thrust::tuple<real_t, real_t> &tpl)
        {
	  const real_t &x = thrust::get<0>(tpl); 
	  const real_t &c = thrust::get<1>(tpl); 

          return x >= c ? y : 0; 
        }
      };
    }  

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_all()
    {
      hskpng_sort(); 

      thrust_device::vector<real_t> &n_filtered(tmp_device_real_part);

      thrust::copy(
        n.begin(), n.end(),
        n_filtered.begin()
      );

      selected_before_counting = true;
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_rng(
      const real_t &min, const real_t &max, 
      const typename thrust_device::vector<real_t>::iterator &vec_bgn
    )
    {
      hskpng_sort(); 

      // transforming n -> n if within range, else 0
      thrust_device::vector<real_t> &n_filtered(tmp_device_real_part);

      thrust::transform(
        n.begin(), n.end(),   // input - 1st arg
	vec_bgn,              // input - 2nd arg
	n_filtered.begin(),   // output
	detail::range_filter<real_t>(min, max) 
      );

      selected_before_counting = true;
    }
 
    // selects particles for which vec1[i] >= vec2[i]
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_cmp(
      const typename thrust_device::vector<real_t>::iterator &vec1_bgn,
      const typename thrust_device::vector<real_t>::iterator &vec2_bgn
    )
    {
      hskpng_sort();

      thrust_device::vector<real_t> &n_filtered(tmp_device_real_part);

      thrust::transform(
        n.begin(), n.end(),                      // input - 1st arg
        thrust::make_zip_iterator(               //
          thrust::make_tuple(vec1_bgn, vec2_bgn) // input - 2nd arg
        ),                                       // 
        n_filtered.begin(),                      // output
        detail::comp_filter<real_t>()            // op
      );

      selected_before_counting = true;
    }

    // selects particles for which vec[i] >= 0
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_ge0(
      const typename thrust_device::vector<real_t>::iterator &vec_bgn
    )
    {
      hskpng_sort();

      thrust_device::vector<real_t> &n_filtered(tmp_device_real_part);
      
      {
        namespace arg = thrust::placeholders;
        thrust::transform(
          n.begin(), n.end(),                      // input - 1st arg
          vec_bgn,                                 // input - 2nd arg
          n_filtered.begin(),                      // output
          arg::_1 * (arg::_2 >= 0)                 // op
        );
      }

      selected_before_counting = true;
    }

    namespace detail
    {
      template <typename real_t>
      struct moment_counter : thrust::unary_function<const thrust::tuple<real_t, real_t>&, real_t>
      {
        real_t xp;

        moment_counter(real_t xp) : xp(xp) {}

        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t> &tpl)
        {
          const real_t n = thrust::get<0>(tpl);
          const real_t x = thrust::get<1>(tpl);
          return n * pow(x, xp); // TODO: check if xp=0 is optimised
        }
      };

      template <typename real_t, typename n_t>
      struct moment_counter_cond : thrust::unary_function<const thrust::tuple<n_t, real_t>&, real_t>
      {
        real_t xp;

        moment_counter_cond(real_t xp) : xp(xp) {}

        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<n_t, real_t> &tpl)
        {
          const real_t n = thrust::get<0>(tpl);
          const real_t x = thrust::get<1>(tpl);
          return n * pow(x, xp); // TODO: check if xp=0 is optimised
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_calc(
      const typename thrust_device::vector<real_t>::iterator &vec_bgn,
      const real_t power
    )
    {
      assert(selected_before_counting);

      // same as above
      thrust_device::vector<real_t> &n_filtered(tmp_device_real_part);

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::const_iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_t;
      typedef thrust::zip_iterator<thrust::tuple<pi_t, pi_t> > zip_it_t;

      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        typename thrust_device::vector<real_t>::iterator
      > n = thrust::reduce_by_key(
        // input - keys
        sorted_ijk.begin(), sorted_ijk.end(),  
        // input - values
        thrust::make_transform_iterator(
	  zip_it_t(thrust::make_tuple(
            pi_t(n_filtered.begin(),   sorted_id.begin()),
            pi_t(vec_bgn,              sorted_id.begin())
          )),
          detail::moment_counter<real_t>(power)
        ),
        // output - keys
        count_ijk.begin(),
        // output - values
        count_mom.begin()
      );  

      count_n = n.first - count_ijk.begin();
      assert(count_n > 0 && count_n <= n_cell);

      // dividing by dv
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,     // input - first arg
        thrust::make_permutation_iterator(                  // input - second arg
          dv.begin(),
	  count_ijk.begin()
        ),
        count_mom.begin(),                                  // output (in place)
        thrust::divides<real_t>()
      );

      // dividing by rhod to get specific moments
      // (for compatibility with blk_1m and blk_2m reporting mixing ratios)
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,     // input - first arg
        thrust::make_permutation_iterator(                  // input - second arg
          rhod.begin(),
	  count_ijk.begin()
        ),
        count_mom.begin(),                                  // output (in place)
        thrust::divides<real_t>()
      );
    }

    // compute non-specific moment for all SDs
    // a separate function since it is needed by condensation
    // TODO: make it a non-specific moment counter for selected SDs,
    // then reuse it in moms_calc (selecting all SDs during condensation
    // doesn't change performance much)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_calc_cond(
      const typename thrust_device::vector<real_t>::iterator &vec_bgn,
      const real_t power
    )
    {
      assert(sorted);

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::const_iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_t;
      typedef thrust::permutation_iterator<
        typename thrust_device::vector<n_t>::const_iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_n_t;
      typedef thrust::zip_iterator<thrust::tuple<pi_n_t, pi_t> > zip_it_t;

      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        typename thrust_device::vector<real_t>::iterator
      > n = thrust::reduce_by_key(
        // input - keys
        sorted_ijk.begin(), sorted_ijk.end(),  
        // input - values
        thrust::make_transform_iterator(
	  zip_it_t(thrust::make_tuple(
            pi_n_t(this->n.begin(),   sorted_id.begin()),
            pi_t(vec_bgn,              sorted_id.begin())
          )),
          detail::moment_counter_cond<real_t, n_t>(power)
        ),
        // output - keys
        count_ijk.begin(),
        // output - values
        count_mom.begin()
      );  

      count_n = n.first - count_ijk.begin();
      assert(count_n > 0 && count_n <= n_cell);
    }
  };  
};
