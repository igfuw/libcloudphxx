// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/count.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      void copy_prop(
        const typename thrust_device::vector<real_t>::iterator &prop_bgn,
        const thrust_device::vector<thrust_size_t> &sorted_id,
        const thrust_size_t &n_flagged
      ) 
      {
	thrust::copy_n(
          // input (properties of n_flagged particles with largest multiplicities
          thrust::make_permutation_iterator(prop_bgn, make_reverse_iterator(sorted_id.end())),
          // count
          n_flagged,
          // output (properties of n_flagged particles with zero multiplicities
          thrust::make_permutation_iterator(prop_bgn, sorted_id.begin())
	);
      }
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::rcyc(
    )
    {   
      // count the numer of paticles to recycle
      thrust_size_t n_flagged;
      {
	namespace arg = thrust::placeholders;
        n_flagged = thrust::count_if(n.begin(), n.end(), arg::_1 == 0);
      }
      assert(n_flagged <= n_part / 2);

      if (n_flagged == 0) return;

      // sort according to multiplicity 
      // -> on one end: those flagged for recycling 
      // -> on the other end: those that will be splitted

      // using sorted_id and sorted_ijk as temporary space - anyhow, after recycling these are not valid anymore!
      sorted = false;
      thrust::sequence(sorted_id.begin(), sorted_id.end()); 
      {
#if defined(__NVCC__) 
        assert(sizeof(thrust_size_t) == sizeof(n_t));
#else
        static_assert(sizeof(thrust_size_t) == sizeof(n_t), "");
#endif
	thrust_device::vector<thrust_size_t> &tmp(sorted_ijk);
	thrust::copy(n.begin(), n.end(), tmp.begin());

	thrust::sort_by_key(
	  tmp.begin(), tmp.end(),
	  sorted_id.begin()
	);
      }

      // for each property... 
      if (opts_init.nx > 0) detail::copy_prop<real_t>(x.begin(), sorted_id, n_flagged); 
      if (opts_init.ny > 0) detail::copy_prop<real_t>(y.begin(), sorted_id, n_flagged); 
      if (opts_init.nz > 0) detail::copy_prop<real_t>(z.begin(), sorted_id, n_flagged); 

      detail::copy_prop<real_t>(rd3.begin(), sorted_id, n_flagged);
      detail::copy_prop<real_t>(rw2.begin(), sorted_id, n_flagged);
      detail::copy_prop<real_t>(kpa.begin(), sorted_id, n_flagged);

      // TODO: only if chem enabled
      for (int i = 0; i < chem_aq_n; ++i)
        detail::copy_prop<real_t>(chem_bgn[i], sorted_id, n_flagged);

      {
        namespace arg = thrust::placeholders;

        // increasing multiplicities of recycled particles
	thrust::transform(
	  // input 
          thrust::make_permutation_iterator(n.begin(), make_reverse_iterator(sorted_id.end())),
          thrust::make_permutation_iterator(n.begin(), make_reverse_iterator(sorted_id.end())) + n_flagged,
	  // output
          thrust::make_permutation_iterator(n.begin(), sorted_id.begin()),
          // op
          arg::_1 - (arg::_1 / 2)
	);

	// reducing multiplicites of splitted particles
	thrust::transform(
          // input
          thrust::make_permutation_iterator(n.begin(), make_reverse_iterator(sorted_id.end())),
          thrust::make_permutation_iterator(n.begin(), make_reverse_iterator(sorted_id.end())) + n_flagged,
          // output
          thrust::make_permutation_iterator(n.begin(), make_reverse_iterator(sorted_id.end())),
          // op
          arg::_1 / 2
	);
      };
    }
  };  
};
