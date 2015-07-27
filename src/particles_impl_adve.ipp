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
      struct adve_helper
      {
        const real_t dx;

        adve_helper(const real_t dx) : dx(dx) {} 

        BOOST_GPU_ENABLED
        real_t operator()(thrust::tuple<real_t, const thrust_size_t, const real_t, const real_t> tpl)
        {
          real_t x = thrust::get<0>(tpl);
          const thrust_size_t floor_x_over_dx = thrust::get<1>(tpl);
          real_t C_l = thrust::get<2>(tpl);
          real_t C_r = thrust::get<3>(tpl);

          // integrating using backward Euler scheme + interpolation/extrapolation
          // 
          // x_new = x_old + v(x_new) * dt = x_old + C(x_new) * dx
          //    
          //     C(x_new) = (1-w) * C_l + w * C_r 
          //     w = x_new/dx - floor(x_old/dx) 
          //
          // x_new = x_old + C_l * dx + w * dx * (C_r - C_l)
          //       = x_old + C_l * dx + x_new * (C_r - C_l) - dx * floor(x_old/dx) * (C_r - C_l)
          // 
          // x_new * (1 - C_r + C_l) = x_old + C_l * dx - dx * floor(x_old/dx) * (C_r - C_l)
          // x_new =  (x_old + C_l * dx - dx * floor(x_old/dx) * (C_r - C_l)) / (1 - C_r + C_l)

          return (
            x + dx * (C_l - floor_x_over_dx * (C_r - C_l))
          ) / (
            1 - (C_r - C_l)
          );
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::adve()
    {   
      switch (n_dims)
      {
        case 3:
        case 2:
        {
          // advection
          typedef thrust::permutation_iterator<
            typename thrust_device::vector<thrust_size_t>::iterator, 
            typename thrust_device::vector<thrust_size_t>::iterator
          > pi_size_size;
          typedef thrust::permutation_iterator<
            typename thrust_device::vector<real_t>::iterator, 
            pi_size_size
          > pi_real_size;

          {
            const pi_real_size
              C_lft(courant_x.begin(), pi_size_size(lft.begin(), ijk.begin())),
              C_rgt(courant_x.begin(), pi_size_size(rgt.begin(), ijk.begin()));
            thrust::transform_if(
              thrust::make_zip_iterator(make_tuple(x.begin(), i.begin(), C_lft,        C_rgt       )), // input - begin
              thrust::make_zip_iterator(make_tuple(x.end(),   i.end(),   C_lft+opts_init.n_sd_max, C_rgt+opts_init.n_sd_max)), // input - end
              sd_stat.begin(), // stencil
              x.begin(), // output
              detail::adve_helper<real_t>(opts_init.dx),
              detail::is_active() // move only active particles
            );
          }

          if (n_dims == 3)
          {
            const pi_real_size
              C_fre(courant_x.begin(), pi_size_size(fre.begin(), ijk.begin())),
              C_hnd(courant_x.begin(), pi_size_size(hnd.begin(), ijk.begin()));
            thrust::transform_if(
              thrust::make_zip_iterator(make_tuple(y.begin(), j.begin(), C_fre,        C_hnd       )), // input - begin
              thrust::make_zip_iterator(make_tuple(y.end(),   j.end(),   C_fre+opts_init.n_sd_max, C_hnd+opts_init.n_sd_max)), // input - end
              sd_stat.begin(), // stencil
              y.begin(), // output
              detail::adve_helper<real_t>(opts_init.dy),
              detail::is_active() // move only active particles
            );
          }

          {
            const pi_real_size
              C_abv(courant_z.begin(), pi_size_size(abv.begin(), ijk.begin())),
              C_blw(courant_z.begin(), pi_size_size(blw.begin(), ijk.begin()));
            thrust::transform_if(
              thrust::make_zip_iterator(make_tuple(z.begin(), k.begin(), C_blw,        C_abv       )), // input - begin
              thrust::make_zip_iterator(make_tuple(z.end(),   k.end(),   C_blw+opts_init.n_sd_max, C_abv+opts_init.n_sd_max)), // input - end
              sd_stat.begin(), // stencil
              z.begin(), // output
              detail::adve_helper<real_t>(opts_init.dz),
              detail::is_active() // move only active particles
            );
          }

          break; 
        }
        case 1:
          assert(false && "TODO");
          break;
        case 0: break;
        default: assert(false);
      }
    }
  };  
};
