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
    void particles_t<real_t, device>::impl::ravel_ijk(const thrust_size_t begin_shift, const bool refined)
    {
      auto & _ijk = refined ? ijk_ref : ijk;
      switch (n_dims)
      {
        case 0: 
          break;
        case 1:
          thrust::copy(i.begin()+begin_shift, i.end(), _ijk.begin()+begin_shift);
          break;
        case 2:
          namespace arg = thrust::placeholders;
          thrust::transform(
            i.begin()+begin_shift, i.end(), // input - first arg
            k.begin()+begin_shift,          // input - second arg
            _ijk.begin()+begin_shift,        // output
            arg::_1 * opts_init.nz + arg::_2   // assuming z varies first
          );
          break;
        case 3:
          namespace arg = thrust::placeholders;
          thrust::transform(
            i.begin()+begin_shift, i.end(), // input - first arg
            j.begin()+begin_shift,          // input - second arg
            _ijk.begin()+begin_shift,        // output
            arg::_1 * (opts_init.nz * opts_init.ny) + 
            arg::_2 * opts_init.nz
          );
          thrust::transform(
            _ijk.begin()+begin_shift, _ijk.end(),
            k.begin()+begin_shift,
            _ijk.begin()+begin_shift, // in-place!
            arg::_1 + arg::_2
          );
          // TODO: replace these two transforms with single one
          break;
        default:
          assert(false);
      }
    }


    // calc i, j and k from ijk
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::unravel_ijk(const thrust_size_t begin_shift) // default = 0
    {
      switch(n_dims)
      {
        case 3:
          namespace arg = thrust::placeholders;
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
        case 2:
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
        case 1:
          thrust::copy(ijk.begin() + begin_shift, ijk.end(), i.begin() + begin_shift); // only x
        case 0:
          break;
        default:
          assert(false);
          break;
      }
    }

    namespace {
      // helper functor to calculate cell index
      template <typename real_t>
      struct cell_idx_hlpr_t{

        const thrust_size_t begin_shift;

        void operator()(
          thrust_device::vector<real_t> &vx,
          thrust_device::vector<thrust_size_t> &vi,
          const real_t &vd
        ) {
          thrust::transform(
            vx.begin() + begin_shift, vx.end(),                                // input
            vi.begin() + begin_shift,                                          // output
            detail::divide_by_constant_and_cast<double, thrust_size_t>(vd) // has to be done on doubles to avoid i==nx due to low precision of nvcc math; TODO: now that rand uniform has range [0,1), float might be good here?
          );
        }

        cell_idx_hlpr_t(const thrust_size_t begin_shift = 0):
          begin_shift(begin_shift)
          {}
      };
    };

    // ijk_ref - cell number in the refined grid
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_ijk_ref(const thrust_size_t begin_shift, const bool _unravel_ijk)
    {   
      if(opts_init.n_ref <= 1) return;

      cell_idx_hlpr_t<real_t> helper(begin_shift);

      if (opts_init.nx != 0) helper(x, i, opts_init.dx / opts_init.n_ref);
      if (opts_init.ny != 0) helper(y, j, opts_init.dy / opts_init.n_ref);
      if (opts_init.nz != 0) helper(z, k, opts_init.dz / opts_init.n_ref);

      // raveling i, j & k into ijk_ref
      ravel_ijk(begin_shift, true);

      // go back to i, j, k in the normal cell range
      if(_unravel_ijk)
        unravel_ijk(begin_shift);
    }   

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_ijk()
    {   
      hskpng_ijk_ref(0, false);

      cell_idx_hlpr_t<real_t> helper;

      if (opts_init.nx != 0) helper(x, i, opts_init.dx);
      if (opts_init.ny != 0) helper(y, j, opts_init.dy);
      if (opts_init.nz != 0) helper(z, k, opts_init.dz);

      // raveling i, j & k into ijk
      ravel_ijk();
      
      // flagging that particles are no longer sorted 
      sorted = false;
    }   
  };  
};
