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
	const opts_init_t<real_t> &o;

        dv_eval(const opts_init_t<real_t> &o) : o(o) {}

        __device__
        real_t operator()(const int &ijk)
        {
#if !defined(__NVCC__)
          using std::min;
          using std::max;
#endif

          // ijk = (i*nj + j)*nk + k
          const int
            i = (ijk / max(1,o.nz)) / max(1,o.ny),
            j = (ijk / max(1,o.nz)) % max(1,o.ny),
            k =  ijk % max(1,o.nz);
             
          return 
	    (min((i + 1) * o.dx, o.x1) - max(i * o.dx, o.x0)) *
	    (min((j + 1) * o.dy, o.y1) - max(j * o.dy, o.y0)) *
	    (min((k + 1) * o.dz, o.z1) - max(k * o.dz, o.z0));
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_grid()
    {
      using namespace thrust::placeholders;

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
            _1
	  );
	  thrust::transform(
            zero, zero + n_cell, // input - 1st arg
            rgt.begin(),         // output
            _1 + opts_init.nx
	  );
	  thrust::transform(
            zero, zero + n_cell, // input - 1st arg
            blw.begin(),         // output
            _1 + (_1 / opts_init.nz)
	  );
	  thrust::transform(
            zero, zero + n_cell, // input - 1st arg
            abv.begin(),         // output
            _1 + (_1 / opts_init.nz) + 1
	  );

	  break;
        case 0: break;
	default: assert(false && "TODO");
      }
    }
  };
};
