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
        throw std::runtime_error("SDs were to be removed, but it is not yet compatible with chemistry");

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
            thrust::make_zip_iterator(thrust::make_tuple(z.begin(), k.begin()))
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
      {
        thrust_device::vector<real_t> *vec[] = {&rw2, &rd3, &kpa, &x, &y, &z, &vt, &tmp_device_real_part, &u01};
        for(int i=0; i<9; ++i)
          vec[i]->erase(vec[i]->begin() + n_part, vec[i]->end());
      }
      {
        thrust_device::vector<thrust_size_t> *vec[] = {&i, &j, &k, &ijk, &sorted_id, &sorted_ijk};
        for(int i=0; i<6; ++i)
          vec[i]->erase(vec[i]->begin() + n_part, vec[i]->end());
      }
      n.erase(n.begin() + n_part, n.end());
      un.erase(un.begin() + n_part, un.end());
    }
  };  
};
