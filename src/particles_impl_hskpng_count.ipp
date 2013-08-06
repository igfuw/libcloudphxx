// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/reduce.h>
#include <thrust/iterator/constant_iterator.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, int device>
    void particles<real_t, device>::impl::hskpng_count()
    {   
      hskpng_sort();

      // computing count_* - number of particles per grid cell
      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        thrust_device::vector<n_t>::iterator
      > n = thrust::reduce_by_key(
        sorted_ijk.begin(), sorted_ijk.end(),   // input - keys
        thrust::make_constant_iterator(n_t(1)), // input - values
        count_ijk.begin(),                      // output - keys
        count_num.begin()                       // output - values
      );
      count_n = n.first - count_ijk.begin();
      assert(count_n > 0 && count_n <= n_cell);
    }   


// TODO: move to another file?
    namespace detail
    {
      template <typename real_t>
      struct moment_counter_dry
      {
        __device__
        thrust::tuple<real_t> operator()()
        {
          return thrust::make_tuple(real_t(666));
        }
      };
    }  

    template <typename real_t, int device>
    void particles<real_t, device>::impl::hskpng_rd_moms()
    {
      hskpng_sort();

/* TODO!!!
      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        thrust_device::vector<n_t>::iterator
      > n = thrust::reduce_by_key(
        // input - keys
        sorted_ijk.begin(), sorted_ijk.end(),  
        // input - values
        thrust::transform_iterator<    
          detail::moment_counter_dry<real_t>,
          thrust_device::vector<thrust_size_t>::iterator,
          thrust::tuple<real_t> // TODO: other moments?
        >(sorted_id.begin(), detail::moment_counter_dry<real_t>()),
        // output - keys
        count_ijk.begin(),
        // output - values
        thrust::make_zip_iterator(thrust::make_tuple(count_num))
      );  

      count_n = n.first - count_ijk.begin();
      assert(count_n > 0 && count_n <= n_cell);
*/
    }
  };  
};
