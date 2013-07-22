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
    void particles<real_t, device>::impl::init_e2l(const ptrdiff_t *strides)
    {
      // memory allocation
      rhod.resize(n_cell);
      rhod_th.resize(n_cell);
      rhod_rv.resize(n_cell);
      l2e.resize(n_cell);
    
      // filling in l2e with values
      switch (n_dims)
      {
	using namespace thrust::placeholders;
	case 0:  
	  l2e[0] = 0;  
	  break;
	case 1:  
	  thrust::transform(zero, zero + n_cell, l2e.begin(), strides[0] * _1); 
	  break;
	case 2:
          // z likely veries fastest as parallelising over x is simpler
	  thrust::transform(zero, zero + n_cell, l2e.begin(), 
	    strides[1] * /* i = */ (_1 % opts.nz) + 
	    strides[0] * /* j = */ (_1 / opts.nz)  
	  );
	  break;
        // TODO: 3D case
	default: assert(false);
      }
    }
  };
};
