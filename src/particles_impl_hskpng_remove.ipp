// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/remove.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
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
    };

    // remove SDs with n=0
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_remove_n0()
    {
      typedef thrust::detail::normal_iterator<thrust_device::pointer<real_t> > it_real_t;
      typedef thrust::detail::normal_iterator<thrust_device::pointer<n_t> > it_n_t;
      typedef thrust::detail::normal_iterator<thrust_device::pointer<thrust_size_t> > it_thrust_size_t;
      typedef thrust::tuple<it_n_t, it_real_t, it_real_t, it_real_t, it_real_t, it_thrust_size_t> tup_params_t;
 
      tup_params_t tup_params = thrust::make_tuple(n.begin(), rw2.begin(), rd3.begin(), kpa.begin(), vt.begin(), ijk.begin());

      if(opts_init.chem_switch)
      {
        // TODO: remove all chem in one remove_if call
        for (int i = chem_all-1; i >= 0; --i)
        {
          namespace arg = thrust::placeholders;

          typename thrust_device::vector<real_t>::iterator new_last = thrust::remove_if(
            chem_bgn[i],
            chem_end[i],
            n.begin(),
            arg::_1 == 0
          );
          
          thrust_device::vector<real_t> &vec(
            i < chem_rhs_beg 
              ? chem_ante_rhs
              : i < chem_rhs_fin
                ? chem_rhs
                : chem_post_rhs
          );
 
          vec.erase(new_last, chem_end[i]);
        }
      }

      if(n_dims == 3)
      {
        typedef thrust::zip_iterator<
          thrust::tuple<
            thrust::zip_iterator<tup_params_t>,
            thrust::zip_iterator<
              thrust::tuple<
                it_real_t, it_real_t, it_real_t,                    //x, y, z
                it_thrust_size_t,it_thrust_size_t,it_thrust_size_t  // i, j, k
              >
            >
          >
        > zip_param_pos_t;
       
        zip_param_pos_t zip_param_pos(
          thrust::make_tuple(
            thrust::make_zip_iterator(tup_params), 
            thrust::make_zip_iterator(thrust::make_tuple(x.begin(), y.begin(), z.begin(), i.begin(), j.begin(), k.begin()))
          )
        );

        zip_param_pos_t new_end = thrust::remove_if(
          zip_param_pos,
          zip_param_pos + n_part, 
          detail::n_eq_zero()
        );
        n_part = new_end - zip_param_pos;
      }
      else if(n_dims == 2)
      {
        typedef thrust::zip_iterator<
          thrust::tuple<
            thrust::zip_iterator<tup_params_t>,
            thrust::zip_iterator<
              thrust::tuple<
                it_real_t, it_real_t,              //x, z
                it_thrust_size_t,it_thrust_size_t  // i, k
              >
            >
          >
        > zip_param_pos_t;
       
        zip_param_pos_t zip_param_pos(
          thrust::make_tuple(
            thrust::make_zip_iterator(tup_params), 
            thrust::make_zip_iterator(thrust::make_tuple(x.begin(), z.begin(), i.begin(), k.begin()))
          )
        );

        zip_param_pos_t new_end = thrust::remove_if(
          zip_param_pos,
          zip_param_pos + n_part, 
          detail::n_eq_zero()
        );
        n_part = new_end - zip_param_pos;
      }
      else if(n_dims == 1)
      {
        typedef thrust::zip_iterator<
          thrust::tuple<
            thrust::zip_iterator<tup_params_t>,
            thrust::zip_iterator<
              thrust::tuple<
                it_real_t,        // z
                it_thrust_size_t  // k
              >
            >
          >
        > zip_param_pos_t;
       
        zip_param_pos_t zip_param_pos(
          thrust::make_tuple(
            thrust::make_zip_iterator(tup_params), 
            thrust::make_zip_iterator(thrust::make_tuple(x.begin(), i.begin()))
          )
        );

        zip_param_pos_t new_end = thrust::remove_if(
          zip_param_pos,
          zip_param_pos + n_part, 
          detail::n_eq_zero()
        );
        n_part = new_end - zip_param_pos;
      }
      else if(n_dims == 0)
      {
        typedef thrust::zip_iterator<tup_params_t> zip_param_pos_t;
       
        zip_param_pos_t zip_param_pos(
            thrust::make_zip_iterator(tup_params)
        );

        zip_param_pos_t new_end = thrust::remove_if(
          zip_param_pos,
          zip_param_pos + n_part, 
          detail::n_eq_zero_0D()
        );
        n_part = new_end - zip_param_pos;
      }

      // resize vectors
      hskpng_resize_npart();

      // resize chem vectors and update chem iterators
      if(opts_init.chem_switch)
        init_chem();
    }
  };  
};
