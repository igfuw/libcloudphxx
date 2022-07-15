// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <thrust/sequence.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class res_t>
      struct arbitrary_sequence //fill container with n 0s, m 1s, l 2s, etc...
      {
        thrust_device::pointer<res_t> res;
        arbitrary_sequence(thrust_device::pointer<res_t> res): res(res) {}
      
        template<typename Tuple>
        BOOST_GPU_ENABLED
        void operator()(Tuple tup)
        {
          for(int i=0; i<thrust::get<0>(tup); ++i)
            *(res+i+thrust::get<1>(tup)) = thrust::get<2>(tup);
        }
      };
    };

    // init ijk based on the number of SDs to init in each cell stored in count_num
    // Particles to init are considered to be sorted by cell number, in order
    // to obtain uniform initial distribution in each cell (see particles_impl_init_dry)
    // reused in source
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_ijk()
    {
      thrust_device::vector<thrust_size_t> &ptr(tmp_device_size_cell);
      thrust::exclusive_scan(count_num.begin(), count_num.end(), ptr.begin()); // number of SDs in cells to init up to (i-1)

      // fill ijk with cell number of each SD
      thrust::for_each(
        thrust::make_zip_iterator(thrust::make_tuple(
          count_num.begin(), ptr.begin(), zero
        )), 
        thrust::make_zip_iterator(thrust::make_tuple(
          count_num.end(), ptr.end(), zero + n_cell
        )), 
        detail::arbitrary_sequence<thrust_size_t>(&(ijk.begin() + n_part_old))
      );
    }
  };
};
