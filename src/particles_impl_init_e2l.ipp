// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      struct periodic_cellno
      {
        int n_x_tot; 
        periodic_cellno(const int &n_x_tot): n_x_tot(n_x_tot) {}
        BOOST_GPU_ENABLED
        int operator()(int x_cell)
        {
          if(x_cell <= n_x_tot && x_cell >= 0) // what if courant_x > 1 and halo of 1 is not enough?
            return x_cell;
          else if (x_cell > n_x_tot)
            return 0;
          else 
            return n_x_tot;
        }
      };
    };
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_e2l(
      const arrinfo_t<real_t> &arr,
      thrust_device::vector<real_t> * key,
      const int ext_x, const int ext_y, const int ext_z,
      const int offset
    )
    {
      // allocating and filling in l2e with values
      l2e[key].resize(key->size());
      switch (n_dims)
      {
	namespace arg = thrust::placeholders;
	case 0:  
	  l2e[key][0] = 0;  
	  break;
	case 1:
        {
          assert(arr.strides[0] == 1);
          int shift =    // index of element of arr copied to 0-th position in key
            + n_cell_bfr // cells in other memory
            + offset;    // additional cells in other memory for arrays bigger than nx*ny*nz like courant numbers
                         // or halo
	  thrust::transform(
            // input
            zero + shift,
            zero + shift + l2e[key].size(), 
            // output
            l2e[key].begin(), 
            // op
            detail::periodic_cellno(n_x_tot) // apply periodic bcnd (if = -1 then = glob_nx; if = glob_nx+1 then = 0)
	  );
	  break;
        }
	case 2:
          // assumes z veries fastest
          assert(arr.strides[1] == 1);
	  thrust::transform(
            // input
            zero + n_cell_bfr + offset, zero + n_cell_bfr + offset + l2e[key].size(), 
            // output
            l2e[key].begin(), 
            // op
	    arr.strides[0] * /* i = */ (arg::_1 / (opts_init.nz + ext_z)) +
	    arr.strides[1] * /* j = */ (arg::_1 % (opts_init.nz + ext_z))  
	  );
	  break;
        case 3:
          assert(arr.strides[2] == 1);
          thrust::transform(
            // input
            zero + n_cell_bfr + offset, zero + n_cell_bfr + offset + l2e[key].size(),
            // output
            l2e[key].begin(),
            // op
	    arr.strides[0] * /* i = */ (arg::_1 / ((opts_init.nz + ext_z) * (opts_init.ny + ext_y))) +  
            arr.strides[1] * /* j = */ ((arg::_1 / (opts_init.nz + ext_z)) % (opts_init.ny + ext_y)) + 
	    arr.strides[2] * /* k = */ (arg::_1 % ((opts_init.nz + ext_z)))    
          );
          break;
	default: assert(false);
      }
    }
  };
};
