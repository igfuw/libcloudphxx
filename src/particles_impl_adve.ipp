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
      struct adve_2d
      {
        const real_t dx;

        adve_2d(const real_t dx) : dx(dx) {} 

        __device__
        real_t operator()(thrust::tuple<real_t, const thrust_size_t, const real_t, const real_t, const real_t> tpl)
        {
          real_t x = thrust::get<0>(tpl);
          const thrust_size_t floor_x_over_dx = thrust::get<1>(tpl);
          real_t C_l = thrust::get<2>(tpl) / thrust::get<4>(tpl);
          real_t C_r = thrust::get<3>(tpl) / thrust::get<4>(tpl);

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

      template <typename real_t>
      struct periodic
      { 
        real_t a, b;
        periodic(real_t a, real_t b) : a(a), b(b) {}
        __device__
        real_t operator()(real_t x)
        {
          return a + fmod((x-a) + (b-a), b-a); // this should call CUDA's fmod!
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::adve()
    {   
      switch (n_dims)
      {
        case 3:
          assert(false && "TODO");
          break; 
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

          const pi_real_size
            C_lft(rhod_courant_x.begin(), pi_size_size(lft.begin(), ijk.begin())),
            C_rgt(rhod_courant_x.begin(), pi_size_size(rgt.begin(), ijk.begin())),
            C_abv(rhod_courant_z.begin(), pi_size_size(abv.begin(), ijk.begin())),
            C_blw(rhod_courant_z.begin(), pi_size_size(blw.begin(), ijk.begin()));

          const thrust::permutation_iterator<
            typename thrust_device::vector<real_t>::iterator,
            typename thrust_device::vector<thrust_size_t>::iterator
          > rhod_ijk(rhod.begin(), ijk.begin());

          thrust::transform(
            thrust::make_zip_iterator(make_tuple(x.begin(), i.begin(), C_lft,        C_rgt       , rhod_ijk       )), // input - begin
            thrust::make_zip_iterator(make_tuple(x.end(),   i.end(),   C_lft+n_part, C_rgt+n_part, rhod_ijk+n_part)), // input - end
            x.begin(), // output
            detail::adve_2d<real_t>(opts_init.dx)
          );
          thrust::transform(
            thrust::make_zip_iterator(make_tuple(z.begin(), k.begin(), C_blw,        C_abv       , rhod_ijk       )), // input - begin
            thrust::make_zip_iterator(make_tuple(z.end(),   k.end(),   C_blw+n_part, C_abv+n_part, rhod_ijk+n_part)), // input - end
            z.begin(), // output
            detail::adve_2d<real_t>(opts_init.dz)
          );

          // hardcoded periodic boundary in x! (TODO - not here?)
          thrust::transform(
            x.begin(), x.end(),
            x.begin(),
            detail::periodic<real_t>(opts_init.x0, opts_init.x1)
          );
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
