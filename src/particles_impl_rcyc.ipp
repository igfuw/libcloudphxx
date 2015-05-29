// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/count.h>
#include <thrust/remove.h>

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
 
      struct n_eq_zero
      {
        template <typename Tuple>
        BOOST_GPU_ENABLED
        bool operator()(Tuple tup)
        {
          return(thrust::get<0>(thrust::get<0>(tup)) == 0);
        }      
      };

      struct n_eq_zero_0D
      {
        template <typename Tuple>
        BOOST_GPU_ENABLED
        bool operator()(Tuple tup)
        {
          return(thrust::get<0>(tup) == 0);
        }      
      };

      // TODO: move it somewhere else
      template<typename type>
      void resize_and_free(thrust_device::vector<type> *vec, const thrust_size_t &n)
      {
        if(vec->size() > n) // to deal with unallocated vectors, e. g. in 1-D case
        {
          vec->resize(n);
//          vec->shrink_to_fit(); // should be used to free memory, but crashes with sorted_id...
        }
      }
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::rcyc()
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
       
        // see if there are any SDs to split, if not - remove SDs with n=0
        // same if const multiplicity option is used TODO: in that case sorting above is not necessary
        if(tmp.back()==1 || opts_init.sd_const_multi > 0)
        {
          typedef thrust::detail::normal_iterator<thrust_device::pointer<real_t> > it_real_t;
          typedef thrust::detail::normal_iterator<thrust_device::pointer<n_t> > it_n_t;
          typedef thrust::detail::normal_iterator<thrust_device::pointer<thrust_size_t> > it_thrust_size_t;
          thrust::tuple<it_n_t, it_real_t, it_real_t, it_real_t, it_real_t, it_thrust_size_t> tup_params = thrust::make_tuple(n.begin(), rw2.begin(), rd3.begin(), kpa.begin(), vt.begin(), ijk.begin());

          if(n_dims == 3)
          {
            thrust::remove_if(
	      thrust::make_zip_iterator(thrust::make_tuple(thrust::make_zip_iterator(tup_params), thrust::make_zip_iterator(thrust::make_tuple(x.begin(), y.begin(), z.begin(), i.begin(), j.begin(), k.begin())))),
	      thrust::make_zip_iterator(thrust::make_tuple(thrust::make_zip_iterator(tup_params), thrust::make_zip_iterator(thrust::make_tuple(x.begin(), y.begin(), z.begin(), i.begin(), j.begin(), k.begin())))) + n_part,
              detail::n_eq_zero()
            );
          }
          else if(n_dims == 2)
          {
            thrust::remove_if(
       	      thrust::make_zip_iterator(thrust::make_tuple(thrust::make_zip_iterator(tup_params), thrust::make_zip_iterator(thrust::make_tuple(x.begin(), z.begin(), i.begin(), k.begin())))),
  	      thrust::make_zip_iterator(thrust::make_tuple(thrust::make_zip_iterator(tup_params), thrust::make_zip_iterator(thrust::make_tuple(x.begin(), z.begin(), i.begin(), k.begin())))) + n_part,
              detail::n_eq_zero()
            );
          }
          else if(n_dims == 1)
          {
            thrust::remove_if(
              thrust::make_zip_iterator(thrust::make_tuple(thrust::make_zip_iterator(tup_params), thrust::make_zip_iterator(thrust::make_tuple(z.begin(), k.begin())))),
              thrust::make_zip_iterator(thrust::make_tuple(thrust::make_zip_iterator(tup_params), thrust::make_zip_iterator(thrust::make_tuple(z.begin(), k.begin())))) + n_part,
              detail::n_eq_zero()
            );
          }
          else if(n_dims == 0)
            thrust::remove_if(
  	      thrust::make_zip_iterator(tup_params),
  	      thrust::make_zip_iterator(tup_params) + n_part,
              detail::n_eq_zero_0D()
            );

          n_part -= n_flagged; 
 
          // resize vectors and free memory
          {
            thrust_device::vector<real_t> *vec[] = {&rw2, &rd3, &kpa, &x, &y, &z, &vt, &tmp_device_real_part};
            for(int i=0; i<8; ++i)
              detail::resize_and_free(vec[i],n_part);
          }
          {
            thrust_device::vector<thrust_size_t> *vec[] = {&i, &j, &k, &ijk, &sorted_id, &sorted_ijk};
            for(int i=0; i<6; ++i)
              detail::resize_and_free(vec[i],n_part);
          }
          detail::resize_and_free(&n,n_part);

          return;
        }
      }

      // for each property... 
      if (opts_init.nx > 0) detail::copy_prop<real_t>(x.begin(), sorted_id, n_flagged); 
      if (opts_init.ny > 0) detail::copy_prop<real_t>(y.begin(), sorted_id, n_flagged); 
      if (opts_init.nz > 0) detail::copy_prop<real_t>(z.begin(), sorted_id, n_flagged); 

      detail::copy_prop<real_t>(rd3.begin(), sorted_id, n_flagged);
      detail::copy_prop<real_t>(rw2.begin(), sorted_id, n_flagged);
      detail::copy_prop<real_t>(kpa.begin(), sorted_id, n_flagged);

      // ... chemical properties only if chem enabled
      if (opts_init.chem_switch){
        for (int i = 0; i < chem_aq_n; ++i)
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
    }
  };  
};
