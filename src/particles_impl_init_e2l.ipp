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
    template <typename real_t, int device>
    void particles<real_t, device>::impl::init_e2l(
      const arrinfo_t<real_t> &arr,
      thrust_device::vector<real_t> * key,
      const int ext_x, const int ext_y, const int ext_z
    )
    {
      // allocating and filling in l2e with values
      l2e[key].resize(key->size());
      switch (n_dims)
      {
	using namespace thrust::placeholders;
	case 0:  
	  l2e[key][0] = 0;  
	  break;
	case 1:  
          assert(arr.strides[0] == 1);
          assert(false && "TODO");
	  break;
	case 2:
          // assumes z veries fastest
          assert(arr.strides[1] == 1);
	  thrust::transform(zero, zero + l2e[key].size(), l2e[key].begin(), 
	    arr.strides[0] * /* i = */ (_1 / (opts.nz+ext_z)) +
	    arr.strides[1] * /* j = */ (_1 % (opts.nz+ext_z))  
	  );
	  break;
        case 3:
          assert(arr.strides[2] == 1);
          assert(false && "TODO");
        // TODO: 3D case
	default: assert(false);
      }
    }
  };
};
