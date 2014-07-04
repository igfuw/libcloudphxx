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
      template <typename real_t>
      struct dv_eval
      {
        // note: having a copy of opts_init here causes CUDA crashes (alignment problems?)
        const int 
          nx, ny, nz;
	const real_t 
          dx, dy, dz,
          x0, y0, z0,
          x1, y1, z1;

        dv_eval(const opts_init_t<real_t> &o) : 
          nx(o.nx), ny(o.ny), nz(o.nz),
          dx(o.dx), dy(o.dy), dz(o.dz),
          x0(o.x0), y0(o.y0), z0(o.z0),
          x1(o.x1), y1(o.y1), z1(o.z1) 
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const int &ijk)
        {
#if !defined(__NVCC__)
          using std::min;
          using std::max;
#endif

          const int
            i = (ijk / max(1,nz)) / max(1,ny),
            j = (ijk / max(1,nz)) % max(1,ny),
            k =  ijk % max(1,nz);
          assert(ijk == (i*max(1,ny) + j)*max(1,nz) + k);
             
          return 
	    (min((i + 1) * dx, x1) - max(i * dx, x0)) *
	    (min((j + 1) * dy, y1) - max(j * dy, y0)) *
	    (min((k + 1) * dz, z1) - max(k * dz, z0));
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_grid()
    {
      namespace arg = thrust::placeholders;

      // filling in sample volume data
      dv.resize(n_cell);

      thrust::transform(
	zero, zero + n_cell, // input - 1st arg
	dv.begin(),          // output  
	detail::dv_eval<real_t>(opts_init)
      );

      switch (n_dims)
      {
	case 2:
          // memory allocation
          lft.resize(n_cell);
          rgt.resize(n_cell);
          abv.resize(n_cell);
          blw.resize(n_cell);

          // filling in neighbour info data
	  thrust::transform(
            zero, zero + n_cell, // input - 1st arg
            lft.begin(),         // output
            arg::_1
	  );
	  thrust::transform(
            zero, zero + n_cell, // input - 1st arg
            rgt.begin(),         // output
            arg::_1 + opts_init.nx
	  );
	  thrust::transform(
            zero, zero + n_cell, // input - 1st arg
            blw.begin(),         // output
            arg::_1 + (arg::_1 / opts_init.nz)
	  );
	  thrust::transform(
            zero, zero + n_cell, // input - 1st arg
            abv.begin(),         // output
            arg::_1 + (arg::_1 / opts_init.nz) + 1
	  );

	  break;
        case 0: break;
	default: assert(false && "TODO");
      }
    }
  };
};
