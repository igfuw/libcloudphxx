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
    thrust_size_t particles_t<real_t, device>::impl::rcyc()
    {   

      // count the numer of paticles to recycle
      thrust_size_t n_flagged, n_to_rcyc;
      {
	namespace arg = thrust::placeholders;
        n_flagged = thrust::count_if(n.begin(), n.end(), arg::_1 == 0);
      }

      if (n_flagged == 0) return 0;
      n_to_rcyc = n_flagged;

      if(pure_const_multi) // remove particles if using const_multi, TODO: what if a mixed run, but const_multi prtcls have higher multiplicity? - they will be split
      {
        hskpng_remove_n0();
        return n_flagged;
      }

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
       
        // check how many SDs are available to split
        thrust_size_t n_splittable;
        n_splittable = 
          thrust::find(make_reverse_iterator(tmp.end()), make_reverse_iterator(tmp.begin()), 1) - 
          make_reverse_iterator(tmp.end());

        //if none are splittable remove SDs with n=0
        if(n_splittable==0)
        {
          hskpng_remove_n0();
          return n_to_rcyc;
        }

        // if there are not enough SDs to split, reduce n_flagged
        if(n_splittable < n_flagged) n_flagged = n_splittable;

      }

      // for each property... 
      if (opts_init.nx > 0) detail::copy_prop<real_t>(x.begin(), sorted_id, n_flagged); 
      if (opts_init.ny > 0) detail::copy_prop<real_t>(y.begin(), sorted_id, n_flagged); 
      if (opts_init.nz > 0) detail::copy_prop<real_t>(z.begin(), sorted_id, n_flagged); 

      detail::copy_prop<real_t>(rd3.begin(), sorted_id, n_flagged);
      detail::copy_prop<real_t>(rw2.begin(), sorted_id, n_flagged);
      detail::copy_prop<real_t>(kpa.begin(), sorted_id, n_flagged);
      if(opts_init.sstp_cond > 1 && opts_init.exact_sstp_cond)
      {
        detail::copy_prop<real_t>(sstp_tmp_rv.begin(), sorted_id, n_flagged);
        detail::copy_prop<real_t>(sstp_tmp_th.begin(), sorted_id, n_flagged);
        detail::copy_prop<real_t>(sstp_tmp_rh.begin(), sorted_id, n_flagged);
        if(const_p)
          detail::copy_prop<real_t>(sstp_tmp_p.begin(), sorted_id, n_flagged);
      }

      // ... chemical properties only if chem enabled
      if (opts_init.chem_switch){
        for (int i = 0; i < chem_all; ++i)
          detail::copy_prop<real_t>(chem_bgn[i], sorted_id, n_flagged);
      }

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
      
      // if not all were recycled, remove those with n==0
      if(n_flagged < n_to_rcyc)  hskpng_remove_n0();
      return n_to_rcyc;
    }
  };  
};
