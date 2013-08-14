// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/vterm.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct common__vterm__vt
      {
        __device__ real_t operator()(
          const real_t &rw2, 
          const thrust::tuple<real_t, real_t> &tpl
       ) {   
#if !defined(__NVCC__)
         using std::sqrt;
#endif
         return common::vterm::vt(
           sqrt(rw2)           * si::metres, // TODO: consider caching rw?
           thrust::get<0>(tpl) * si::kelvins,
           thrust::get<1>(tpl) * si::kilograms / si::cubic_metres
         ) / si::metres_per_second;
       }   
      }; 
    };


    template <typename real_t, int device>
    void particles<real_t, device>::impl::hskpng_vterm_invalid()
    {   
      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_t;
      typedef thrust::zip_iterator<thrust::tuple<pi_t, pi_t> > zip_it_t;

      using namespace thrust::placeholders;

      thrust::transform_if(
        rw2.begin(), rw2.end(),                                 // input - 1st arg
	zip_it_t(thrust::make_tuple(
          thrust::make_permutation_iterator(T.begin(),    ijk.begin()),
          thrust::make_permutation_iterator(rhod.begin(), ijk.begin())
        )),                                                     // input - 2nd arg   
        vt.begin(),                                             // condition argument
	vt.begin(),                                             // output
	detail::common__vterm__vt<real_t>(),
        _1 == real_t(detail::invalid)
      );
    }

    template <typename real_t, int device>
    void particles<real_t, device>::impl::hskpng_vterm_all()
    {   
      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_t;
      typedef thrust::zip_iterator<thrust::tuple<pi_t, pi_t> > zip_it_t;

      thrust::transform(
        rw2.begin(), rw2.end(),                                 // input - 1st arg
	zip_it_t(thrust::make_tuple(
          thrust::make_permutation_iterator(T.begin(),    ijk.begin()),
          thrust::make_permutation_iterator(rhod.begin(), ijk.begin())
        )),                                                     // input - 2nd arg
	vt.begin(),                                             // output
	detail::common__vterm__vt<real_t>()
      );
    }
  };  
};
