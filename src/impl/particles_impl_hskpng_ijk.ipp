// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/theta_dry.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      /// @brief returns ret_t(x/c) 
      template <typename arg_t, typename ret_t>
      struct divide_by_constant_and_cast
      {
        arg_t c;
        divide_by_constant_and_cast(arg_t c) : c(c) {}

        BOOST_GPU_ENABLED
        ret_t operator()(arg_t x) 
        { 
          return ret_t(x/c); 
        }
      };
    };

    // calc ijk from i, j and k
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::ravel_ijk(const thrust_size_t begin_shift) // default = 0
    {
      namespace arg = thrust::placeholders;
      switch (n_dims)
      {
        case 0: 
          break;
        case 1:
        {
          thrust_device::vector<thrust_size_t> &i(i_gp->get());
          thrust::copy(i.begin()+begin_shift, i.end(), ijk.begin()+begin_shift);
          break;
        }
        case 2:
        {
          thrust_device::vector<thrust_size_t> &i(i_gp->get());
          thrust_device::vector<thrust_size_t> &k(k_gp->get());
          thrust::transform(
            i.begin()+begin_shift, i.end(), // input - first arg
            k.begin()+begin_shift,          // input - second arg
            ijk.begin()+begin_shift,        // output
            arg::_1 * opts_init.nz + arg::_2   // assuming z varies first
          );
          break;
        }
        case 3:
        {
          thrust_device::vector<thrust_size_t> &i(i_gp->get());
          thrust_device::vector<thrust_size_t> &j(j_gp->get());
          thrust_device::vector<thrust_size_t> &k(k_gp->get());
          thrust::transform(
            i.begin()+begin_shift, i.end(), // input - first arg
            j.begin()+begin_shift,          // input - second arg
            ijk.begin()+begin_shift,        // output
            arg::_1 * (opts_init.nz * opts_init.ny) + 
            arg::_2 * opts_init.nz
          );
          thrust::transform(
            ijk.begin()+begin_shift, ijk.end(),
            k.begin()+begin_shift,
            ijk.begin()+begin_shift, // in-place!
            arg::_1 + arg::_2
          );
          // TODO: replace these two transforms with single one
          break;
        }
        default:
          assert(false);
      }
    }


    // calc i, j and k from ijk
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::unravel_ijk(const thrust_size_t begin_shift) // default = 0
    {
      namespace arg = thrust::placeholders;
      switch(n_dims)
      {
        case 3:
        {
          reset_guardp(i_gp, tmp_device_size_part);  // acquire tmp array to store i
          reset_guardp(j_gp, tmp_device_size_part);  
          reset_guardp(k_gp, tmp_device_size_part);  

          thrust_device::vector<thrust_size_t> &i(i_gp->get());
          thrust_device::vector<thrust_size_t> &j(j_gp->get());
          thrust_device::vector<thrust_size_t> &k(k_gp->get());
          // y
          thrust::transform(
            ijk.begin() + begin_shift, ijk.end(), // input - first arg
            j.begin() + begin_shift,        // output
            (arg::_1 / opts_init.nz) % (opts_init.ny) // z varies first
          );
          // z
          thrust::transform(
            ijk.begin() + begin_shift, ijk.end(), // input - first arg
            k.begin() + begin_shift,        // output
            arg::_1 % (opts_init.nz)   // z varies first
          );
          // x
          thrust::transform(
            ijk.begin() + begin_shift, ijk.end(), // input - first arg
            i.begin() + begin_shift,        // output
            arg::_1 / (opts_init.nz * opts_init.ny)    // z and y vary first
          );
          break;
        }
        case 2:
        {
          reset_guardp(i_gp, tmp_device_size_part);  // acquire tmp array to store i
          reset_guardp(k_gp, tmp_device_size_part);  

          thrust_device::vector<thrust_size_t> &i(i_gp->get());
          thrust_device::vector<thrust_size_t> &k(k_gp->get());
          // z
          thrust::transform(
            ijk.begin() + begin_shift, ijk.end(), // input - first arg
            k.begin() + begin_shift,        // output
            arg::_1 % (opts_init.nz)   // z varies first
          );
          // x
          thrust::transform(
            ijk.begin() + begin_shift, ijk.end(), // input - first arg
            i.begin() + begin_shift,        // output
            arg::_1 / (opts_init.nz)
          );
          break;
        }
        case 1:
        {
          reset_guardp(i_gp, tmp_device_size_part);  // acquire tmp array to store i
          thrust_device::vector<thrust_size_t> &i(i_gp->get());

          thrust::copy(ijk.begin() + begin_shift, ijk.end(), i.begin() + begin_shift); // only x
          break;
        }
        case 0:
          break;
        default:
          assert(false);
          break;
      }
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_ijk()
    {   
      // helper functor
      struct {
        void operator()(
          thrust_device::vector<real_t> &vx,
          thrust_device::vector<thrust_size_t> &vi,
          const real_t &vd
        ) {
          thrust::transform(
            vx.begin(), vx.end(),                                // input
            vi.begin(),                                          // output
            detail::divide_by_constant_and_cast<double, thrust_size_t>(vd) // has to be done on doubles to avoid i==nx due to low precision of nvcc math; TODO: now that rand uniform has range [0,1), float might be good here?
          );
        }
      } helper;

      if (opts_init.nx != 0) 
      {
        reset_guardp(i_gp, tmp_device_size_part);  // acquire tmp array to store i
        thrust_device::vector<thrust_size_t> &i(i_gp->get());
        helper(x, i, opts_init.dx);
      }
      if (opts_init.ny != 0) 
      {
        reset_guardp(j_gp, tmp_device_size_part);  // acquire tmp array to store j
        thrust_device::vector<thrust_size_t> &j(j_gp->get());
        helper(y, j, opts_init.dy);
      }
      if (opts_init.nz != 0) 
      {
        reset_guardp(k_gp, tmp_device_size_part);  // acquire tmp array to store k
        thrust_device::vector<thrust_size_t> &k(k_gp->get());
        helper(z, k, opts_init.dz);
      }

      // raveling i, j & k into ijk
      ravel_ijk();
      
      // flagging that particles are no longer sorted 
      sorted = false;
    }   
  };  
};
