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
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_count()
    {   
      hskpng_sort();


      // computing count_* - number of particles per grid cell
      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        thrust_device::vector<n_t>::iterator
      > it_pair = thrust::reduce_by_key(
        sorted_ijk.begin(), sorted_ijk.end(),   // input - keys
        thrust::make_constant_iterator(n_t(1)), // input - values
        count_ijk.begin(),                      // output - keys
        count_num.begin()                       // output - values
      );
      count_n = it_pair.first - count_ijk.begin();
#if !defined(NDEBUG)
      if(count_n > n_cell)
      {
        std::cerr << "count_n: " << count_n << std::endl;
        std::cerr << "n_cell: " << n_cell << std::endl;
        /*
        std::cerr << "ijk:" << std::endl;
        debug::print(ijk);
        */
        std::cerr << "sorted_ijk:" << std::endl;
        debug::print(sorted_ijk);
        std::cerr << "count_ijk:" << std::endl;
        debug::print(count_ijk);
      }
#endif
      assert(count_n <= n_cell);
    }   
  };  
};
