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
        const long int n_tot;                     // size of the input array
        periodic_cellno(const long int &n_tot):
          n_tot(n_tot) {}
        BOOST_GPU_ENABLED
        long int operator()(long int cell_idx)
        {
          if (cell_idx >= n_tot)
            cell_idx -= n_tot;
          else if (cell_idx < 0)
            cell_idx += n_tot;

          return cell_idx;
        }
      };
    };
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_e2l(
      const arrinfo_t<real_t> &arr,
      thrust::device_vector<real_t> * key,
      const int ext_x, const int ext_y, const int ext_z,
      const long int offset
    )
    {
      // allocating and filling in l2e with values
      l2e[key].resize(key->size());

      long int shift =    // index of element of arr copied to 0-th position in key
        + n_cell_bfr // cells in other memory
        + offset;    // additional cells in other memory for arrays bigger than nx*ny*nz (like courant numbers),
                     // or halo

      switch (n_dims)
      {
	namespace arg = thrust::placeholders;
	case 0:  
	  l2e[key][0] = 0;  
	  break;
	case 1:
          assert(arr.strides[0] == 1);
	  thrust::transform(
            // input
            thrust::make_counting_iterator<int>(0) + shift,                   // long int didnt work
            thrust::make_counting_iterator<int>(0) + shift + l2e[key].size(), 
            // output
            l2e[key].begin(), 
            // op
            arg::_1
	  );
	  break;
	case 2:
          // assume z changes first
          assert(arr.strides[1] == 1);
	  thrust::transform(
            // input
            thrust::make_counting_iterator<int>(0) + shift,
            thrust::make_counting_iterator<int>(0) + shift + l2e[key].size(), 
            // output
            l2e[key].begin(), 
            // op
	    arr.strides[0] * /* i = */ (arg::_1 / (opts_init.nz + ext_z)) +
	    arr.strides[1] * /* j = */ (arg::_1 % (opts_init.nz + ext_z))     // module of negative value might not work in 2003 standard?
	  );
	  break;
        case 3:
          assert(arr.strides[2] == 1);
          thrust::transform(
            // input
            thrust::make_counting_iterator<int>(0) + shift,
            thrust::make_counting_iterator<int>(0) + shift + l2e[key].size(), 
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

      // apply bcnd for halo
      thrust::transform(
        l2e[key].begin(), l2e[key].begin() + l2e[key].size(),
        l2e[key].begin(), // in place 
        detail::periodic_cellno((n_x_tot + ext_x) * arr.strides[0])
      );
    }
  };
};
