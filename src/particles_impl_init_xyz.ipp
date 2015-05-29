// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <thrust/sequence.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      struct arbitrary_sequence //fill container with n 0s, m 1s, l 2s, etc...
      {
        thrust_device::pointer<thrust_size_t> res;
        arbitrary_sequence(thrust_device::pointer<thrust_size_t> res): res(res) {}
      
        template<typename Tuple>
        BOOST_GPU_ENABLED
        void operator()(Tuple tup)
        {
          for(int i=0; i<thrust::get<0>(tup); ++i)
            *(res+i+thrust::get<1>(tup)) = thrust::get<2>(tup);
        }
      };

    };
    // init_xyz, to get uniform distribution in each cell
    // first n_cell SDs are distributed one per each cell,
    // then same with second n_cell particles, etc.
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_xyz()
    {
      thrust_device::vector<real_t> 
                  *v[3] = { &x,           &y,           &z           };
      const int    n[3] = { opts_init.nx, opts_init.ny, opts_init.nz };
      const real_t a[3] = { opts_init.x0, opts_init.y0, opts_init.z0 };
      const real_t d[3] = { opts_init.dx, opts_init.dy, opts_init.dz };
      thrust_device::vector<thrust_size_t> 
                  *ii[3] = { &i,           &j,           &k           };

      if(opts_init.sd_conc > 0)
        thrust::fill(count_num.begin(), count_num.end(), opts_init.sd_conc); // if using const_multi, count_num is already filled

debug::print(count_num);
      thrust_device::vector<thrust_size_t> &ptr(tmp_device_size_cell);
debug::print(ptr);
      thrust::exclusive_scan(count_num.begin(), count_num.end(), ptr.begin()); // number of SDs in cells up to (i-1)
debug::print(count_num);
debug::print(ptr);

      // fill sorted ijk with cell number of each SD
      thrust::for_each(
        thrust::make_zip_iterator(thrust::make_tuple(
          count_num.begin(), ptr.begin(), thrust::make_counting_iterator(0)
        )), 
        thrust::make_zip_iterator(thrust::make_tuple(
          count_num.end(), ptr.end(), thrust::make_counting_iterator(n_cell)
        )), 
        detail::arbitrary_sequence(&(sorted_ijk[0]))
      );
debug::print(sorted_ijk);

      // get i, j, k from sorted_ijk 
      // TODO: check if it is done the same way as syncing with rhod!!
      // i, j, k will be temporarily used
      switch(n_dims)
      {
        case(0):
          break;
        case(1):
          // z
          thrust::copy(sorted_ijk.begin(), sorted_ijk.end(), k.begin());
          break;
        case 2:
          namespace arg = thrust::placeholders;
          // x
          thrust::transform(
            sorted_ijk.begin(), sorted_ijk.end(), // input - first arg
            i.begin(),        // output
            arg::_1 / opts_init.nz   // assuming z varies first
          );
          // z
          thrust::transform(
            sorted_ijk.begin(), sorted_ijk.end(), // input - first arg
            k.begin(),        // output
            arg::_1 % opts_init.nz   // assuming z varies first
          );
          break;
        case 3:
          namespace arg = thrust::placeholders;
          // y
          thrust::transform(
            sorted_ijk.begin(), sorted_ijk.end(), // input - first arg
            j.begin(),        // output
            arg::_1 / (opts_init.nz * opts_init.nx)   // assuming z and x vary first
          );
          // x
          thrust::transform(
            sorted_ijk.begin(), sorted_ijk.end(), // input - first arg
            i.begin(),        // output
            arg::_1 % (opts_init.nz * opts_init.nx) / (opts_init.nz)   // assuming z varies first
          );
          // z
          thrust::transform(
            sorted_ijk.begin(), sorted_ijk.end(), // input - first arg
            k.begin(),        // output
            arg::_1 % (opts_init.nz * opts_init.nx) % (opts_init.nz)   // assuming z varies first
          );
          break;
        default:
          assert(false);
      }
      for (int ix = 0; ix < 3; ++ix)
      {
        if (n[ix] == 0) continue;

        // memory allocation
        v[ix]->resize(n_part);

        // tossing random numbers [0,1] 
        rand_u01(n_part);

	// shifting from [0,1] to random position within respective cell 
        {
          namespace arg = thrust::placeholders;
	  thrust::transform(
	    u01.begin(), 
	    u01.end(),
            ii[ix]->begin(), 
	    v[ix]->begin(), 
	    a[ix] + (arg::_1 + arg::_2) *d[ix]
	  );
        }
      }
    }
  };
};
