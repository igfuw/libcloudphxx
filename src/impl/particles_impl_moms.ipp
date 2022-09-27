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
      const typename thrust_device::vector<real_t>::iterator &vec_bgn,
      const thrust_size_t npart,
      const bool cons // is it a consecutive selection after previous one
    )
    {
      hskpng_sort(); 

      // transforming n -> n if within range, else 0
      thrust_device::vector<real_t> &n_filtered(tmp_device_real_part);

      if(!cons)
        thrust::transform(
          n.begin(), n.begin() + npart, // input - 1st arg
          vec_bgn,                      // input - 2nd arg
          n_filtered.begin(),           // output
          detail::range_filter<real_t>(min, max) 
        );
      else
        thrust::transform(
          n_filtered.begin(), n_filtered.begin() + npart,  // input - 1st arg
          vec_bgn,                                         // input - 2nd arg
          n_filtered.begin(),                              // output
          detail::range_filter<real_t>(min, max) 
        );

      selected_before_counting = true;
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_rng(
      const real_t &min, const real_t &max, 
      const typename thrust_device::vector<real_t>::iterator &vec_bgn,
      const bool cons // is it a consecutive selection after previous one
    )
    {
      moms_rng(min, max, vec_bgn, n_part, cons);
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
#if !defined(__NVCC__)
          using std::pow;
#endif
          const real_t n = thrust::get<0>(tpl);
          const real_t x = thrust::get<1>(tpl);
#if !defined(NDEBUG)
          real_t res;
          if(x >= 0)
            res = n * pow(x, xp);
          else
          {
            int ixp = xp;
            assert(ixp == xp);
            res = n * pow(x, ixp);
          }
          if(isnaninf()(res))
          {
            printf("nan/inf res in moment counter, n = %g x = %g res = %g xp = %g\n", n, x, res, xp);
          }
          return res;
#else
          return x >= 0 ? n * pow(x, xp) : n * pow(x, int(xp)); // for negative x and non-integer xp, pow = NaN
#endif
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_calc(
      const typename thrust_device::vector<real_t>::iterator &vec_bgn,
      const thrust_size_t npart,
      const real_t power,
      const bool specific,
      const bool refined
    )
    {
      assert(selected_before_counting);

      // same as above
      thrust_device::vector<real_t> &n_filtered(tmp_device_real_part);

      thrust_device::vector<thrust_size_t> 
        &_sorted_ijk ( refined ? sorted_ijk.get_ref() : sorted_ijk.get()),
        &_count_ijk  ( refined ? count_ijk.get_ref()  : count_ijk.get());
      thrust_device::vector<real_t> 
        &_count_mom  ( refined ? count_mom.get_ref()  : count_mom.get()),
        &_rhod       ( refined ? rhod.get_ref()       : rhod.get()),
        &_dv         ( refined ? dv.get_ref()         : dv.get());
      auto 
        &_count_n    ( refined ? count_n.get_ref()    : count_n.get()),
        &_n_cell     ( refined ? n_cell.get_ref()     : n_cell.get());

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::const_iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_t;
      typedef thrust::zip_iterator<thrust::tuple<pi_t, pi_t> > zip_it_t;

      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        typename thrust_device::vector<real_t>::iterator
      > it_pair = thrust::reduce_by_key(
        // input - keys
        _sorted_ijk.begin(), _sorted_ijk.begin()+npart,  
        // input - values
        thrust::make_transform_iterator(
          zip_it_t(thrust::make_tuple(
            pi_t(n_filtered.begin(),   sorted_id.begin()),
            pi_t(vec_bgn,              sorted_id.begin())
          )),
          detail::moment_counter<real_t>(power)
        ),
        // output - keys
        _count_ijk.begin(),
        // output - values
        _count_mom.begin()
      );  

      _count_n = it_pair.first - _count_ijk.begin();
#if !defined(NDEBUG)
      {
        int nan_count = thrust::transform_reduce(_count_mom.begin(), _count_mom.begin() + _count_n, isnaninf(), 0, thrust::plus<bool>());
        if(nan_count>0)
        {
          std::cout << nan_count << " nan/inf numbers detected in count_mom after reduce_by_key " << std::endl;
        }
      }
#endif
      assert(_count_n <= _n_cell);
      if(specific)
      {
        // dividing by dv
        thrust::transform(
          _count_mom.begin(), _count_mom.begin() + _count_n,  // input - first arg
          thrust::make_permutation_iterator(                  // input - second arg
            _dv.begin(),
            _count_ijk.begin()
          ),
          _count_mom.begin(),                                  // output (in place)
          thrust::divides<real_t>()
        );
#if !defined(NDEBUG)
        {
          int nan_count = thrust::transform_reduce(_count_mom.begin(), _count_mom.begin() + _count_n, isnaninf(), 0, thrust::plus<bool>());
          if(nan_count>0)
          {
            std::cout << nan_count << " nan/inf numbers detected in count_mom after dividing by dv " << std::endl;
          }
        }
#endif
        // dividing by rhod to get specific moments
        // (for compatibility with blk_1m and blk_2m reporting mixing ratios)
        thrust::transform(
          _count_mom.begin(), _count_mom.begin() + _count_n,     // input - first arg
          thrust::make_permutation_iterator(                     // input - second arg
            _rhod.begin(),                                       
            _count_ijk.begin()
          ),
          _count_mom.begin(),                                  // output (in place)
          thrust::divides<real_t>()
        );
#if !defined(NDEBUG)
        {
          int nan_count = thrust::transform_reduce(_count_mom.begin(), _count_mom.begin() + _count_n, isnaninf(), 0, thrust::plus<bool>());
          if(nan_count>0)
          {
            std::cout << nan_count << " nan/inf numbers detected in count_mom after dividing by rhod " << std::endl;
          }
        }
#endif
      }
#if !defined(NDEBUG)
      {
        int nan_count = thrust::transform_reduce(_count_mom.begin(), _count_mom.begin() + _count_n, isnaninf(), 0, thrust::plus<bool>());
        if(nan_count>0)
        {
          std::cout << nan_count << " nan/inf numbers detected in count_mom " << std::endl;
          std::cout << "count_n:" << _count_n << std::endl;
          std::cout << "count_mom:" << std::endl;
          debug::print(_count_mom.begin(), _count_mom.begin() + _count_n);
          std::cout << "count_ijk:" << std::endl;
          debug::print(_count_ijk.begin(), _count_ijk.begin() + _count_n);
          std::cout << "n_filtered:" << std::endl;
          debug::print(n_filtered);
          std::cout << "sorted_ijk:" << std::endl;
          debug::print(_sorted_ijk);
          std::cout << "sorted_id:" << std::endl;
          debug::print(sorted_id);
          std::cout << "vec:" << std::endl;
          debug::print(vec_bgn, vec_bgn + npart);
          std::cout << "dv:" << std::endl;
          debug::print(_dv);
          std::cout << "rhod:" << std::endl;
          debug::print(_rhod);
        }
      }
#endif
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::moms_calc(
      const typename thrust_device::vector<real_t>::iterator &vec_bgn,
      const real_t power,
      const bool specific,
      const bool refined
    )
    {
      moms_calc(vec_bgn, n_part, power, specific, refined);
    }
  };  
};
