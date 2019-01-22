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

      thrust::device_vector<real_t> &n_filtered(tmp_device_real_part);

      thrust::copy(
        n.begin(), n.end(),
        n_filtered.begin()
      );

      selected_before_counting = true;
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_rng(
      const real_t &min, const real_t &max, 
      const typename thrust::device_vector<real_t>::iterator &vec_bgn
    )
    {
      hskpng_sort(); 

      // transforming n -> n if within range, else 0
      thrust::device_vector<real_t> &n_filtered(tmp_device_real_part);

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
      const typename thrust::device_vector<real_t>::iterator &vec1_bgn,
      const typename thrust::device_vector<real_t>::iterator &vec2_bgn
    )
    {
      hskpng_sort();

      thrust::device_vector<real_t> &n_filtered(tmp_device_real_part);

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
      const typename thrust::device_vector<real_t>::iterator &vec_bgn
    )
    {
      hskpng_sort();

      thrust::device_vector<real_t> &n_filtered(tmp_device_real_part);
      
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
#if !defined(NDEBUG)
          real_t res = n * pow(x, xp); // TODO: check if xp=0 is optimised
          if(isnaninf()(res))
          {
            printf("nan/inf res in moment counter, n = %g x = %g res = %g xp = %g\n", n, x, res, xp);
          }
          return res;
#else
          return n * pow(x, xp); // TODO: check if xp=0 is optimised
#endif
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_calc(
      const typename thrust::device_vector<real_t>::iterator &vec_bgn,
      const real_t power,
      const bool specific
    )
    {
      assert(selected_before_counting);

      // same as above
      thrust::device_vector<real_t> &n_filtered(tmp_device_real_part);

      typedef thrust::permutation_iterator<
        typename thrust::device_vector<real_t>::const_iterator,
        typename thrust::device_vector<thrust_size_t>::iterator
      > pi_t;
      typedef thrust::zip_iterator<thrust::tuple<pi_t, pi_t> > zip_it_t;

      thrust::pair<
        thrust::device_vector<thrust_size_t>::iterator,
        typename thrust::device_vector<real_t>::iterator
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
#if !defined(NDEBUG)
      {
        int nan_count = thrust::transform_reduce(count_mom.begin(), count_mom.begin() + count_n, isnaninf(), 0, thrust::plus<bool>());
        if(nan_count>0)
        {
          std::cout << nan_count << " nan/inf numbers detected in count_mom after reduce_by_key " << std::endl;
        }
      }
#endif
      assert(count_n > 0 && count_n <= n_cell);
      if(specific)
      {
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
#if !defined(NDEBUG)
        {
          int nan_count = thrust::transform_reduce(count_mom.begin(), count_mom.begin() + count_n, isnaninf(), 0, thrust::plus<bool>());
          if(nan_count>0)
          {
            std::cout << nan_count << " nan/inf numbers detected in count_mom after dividing by dv " << std::endl;
          }
        }
#endif
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
#if !defined(NDEBUG)
        {
          int nan_count = thrust::transform_reduce(count_mom.begin(), count_mom.begin() + count_n, isnaninf(), 0, thrust::plus<bool>());
          if(nan_count>0)
          {
            std::cout << nan_count << " nan/inf numbers detected in count_mom after dividing by rhod " << std::endl;
          }
        }
#endif
      }
#if !defined(NDEBUG)
      {
        int nan_count = thrust::transform_reduce(count_mom.begin(), count_mom.begin() + count_n, isnaninf(), 0, thrust::plus<bool>());
        if(nan_count>0)
        {
          std::cout << nan_count << " nan/inf numbers detected in count_mom " << std::endl;
          std::cout << "count_n:" << count_n << std::endl;
          std::cout << "count_mom:" << std::endl;
          debug::print(count_mom.begin(), count_mom.begin() + count_n);
          std::cout << "count_ijk:" << std::endl;
          debug::print(count_ijk.begin(), count_ijk.begin() + count_n);
          std::cout << "n_filtered:" << std::endl;
          debug::print(n_filtered);
          std::cout << "sorted_ijk:" << std::endl;
          debug::print(sorted_ijk);
          std::cout << "sorted_id:" << std::endl;
          debug::print(sorted_id);
          std::cout << "vec:" << std::endl;
          debug::print(vec_bgn, vec_bgn + n_part);
          std::cout << "dv:" << std::endl;
          debug::print(dv);
          std::cout << "rhod:" << std::endl;
          debug::print(rhod);
        }
      }
#endif
    }
  };  
};
