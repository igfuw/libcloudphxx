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
      template <typename real_t, int d>
      struct adve_2d
      {
        thrust_device::vector<real_t> &x;
        const thrust_device::vector<real_t> &C;
        const thrust_device::vector<thrust_size_t> &i, &j, &floor_x;
        const thrust_size_t nz;
        const real_t dx, nxdx;
	const int di, dj;

        adve_2d(
	  thrust_device::vector<real_t> &x,
	  const thrust_device::vector<real_t> &C,
          const real_t dx,
	  const thrust_device::vector<thrust_size_t> &i,
	  const thrust_device::vector<thrust_size_t> &j,
	  const thrust_size_t nx,
	  const thrust_size_t nz
        ) : 
          // copies of the ctor arguments
          C(C), x(x), dx(dx), 
          i(i), j(j), nxdx(nx * dx), nz(nz), 
          // helper vars to work in both dimensions with the same functor
	  di(d == 0 ? 1 : 0),
	  dj(d == 0 ? 0 : 1),
          floor_x(d==0 ? i : j) 
        {} 

        BOOST_GPU_ENABLED
        void operator()(thrust_size_t ix)
        {
          // integrating using backward Euler scheme + interpolation/extrapolation
          // 
          // x_new = x_old + v(x_new) * dt = x_old + C(x_new) * dx
          //    
          //     C(x_new) = (1-w) * C_l + w * C_r 
          //     w = x_new/dx - floor(x_old) 
          //
          // x_new = x_old + C_l * dx + w * dx * (C_r - C_l)
          //       = x_old + C_l * dx + x_new * (C_r - C_l) - floor(x_old) * (C_r - C_l)
          // 
          // x_new * (1 - C_r + C_l) = x_old + C_l * dx - floor(x_old) * (C_r - C_l)
          // x_new =  (x_old + C_l * dx - floor(x_old) * (C_r - C_l)) / (1 - C_r + C_l)
          const real_t 
            C_l = C[(i[ix]   ) * (nz+dj) + j[ix]   ],
            C_r = C[(i[ix]+di) * (nz+dj) + j[ix]+dj];

          x[ix] = (
            x[ix] + C_l * dx - floor_x[ix] * (C_r - C_l)
          ) / (
            1 - C_r + C_l
          );

          // TODO: warning: hardcoded periodic boundary in x!
          if (d==0) x[ix] = fmod(nxdx + x[ix], nxdx);
          
          // TODO: perhaps update i,j,k here?...
        }
      };
    };

    template <typename real_t, int device>
    void particles<real_t, device>::impl::adve()
    {   
      switch (n_dims)
      {
        case 3:
          assert(false && "TODO");
          break; 
        case 2:
	  thrust::for_each(zero, zero + n_part, detail::adve_2d<real_t, 0>(x, courant_x, opts.dx, i, k, opts.nx, opts.nz));
	  thrust::for_each(zero, zero + n_part, detail::adve_2d<real_t, 1>(z, courant_z, opts.dz, i, k, opts.nx, opts.nz));
          break; 
        case 1:
          assert(false && "TODO");
          break;
        case 0: break;
        default: assert(false);
      }
    }
  };  
};
