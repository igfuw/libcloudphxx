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
/*
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
*/
      template <typename real_t>
      struct adve_helper
      {
        const real_t dx;
        bool apply;

        adve_helper(const real_t &dx, bool apply) : dx(dx), apply(apply) {} 

        BOOST_GPU_ENABLED
        real_t operator()(thrust::tuple<real_t, const thrust_size_t, const real_t, const real_t> tpl)
        {
          real_t x = thrust::get<0>(tpl);
          const thrust_size_t floor_x_over_dx = thrust::get<1>(tpl);
          real_t C_l = thrust::get<2>(tpl);
          real_t C_r = thrust::get<3>(tpl);

          // integrating using forward Euler scheme + interpolation/extrapolation
          // 
          // x_new = x_old + v(x_old) * dt = x_old + C(x_old) * dx
          //    
          //     C(x_old) = (1-w) * C_l + w * C_r 
          //     w = x_old/dx - floor(x_old/dx) 
          //
          // x_new = x_old + C_l * dx + w * dx * (C_r - C_l)
          //       = x_old + C_l * dx + x_old * (C_r - C_l) - dx * floor(x_old/dx) * (C_r - C_l)

          // return new position or change in position
          return apply * x + (C_r - C_l) * (x - dx * floor_x_over_dx) + dx * C_l; 
        }
      };
    };
    // calculate change of position due to advection
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::adve_calc(bool apply) // true - save new position, false - save change in position
    {   
      switch (n_dims)
      {
        case 3:
        case 2:
        case 1:
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
            thrust::transform(
              thrust::make_zip_iterator(make_tuple(x.begin(), i.begin(), C_lft,        C_rgt       )), // input - begin
              thrust::make_zip_iterator(make_tuple(x.end(),   i.end(),   C_lft+n_part, C_rgt+n_part)), // input - end
              x.begin(), // output
              detail::adve_helper<real_t>(opts_init.dx, apply)
            );
          }

          if (n_dims > 2)
          {
            const pi_real_size
              C_fre(courant_y.begin(), pi_size_size(fre.begin(), ijk.begin())),
              C_hnd(courant_y.begin(), pi_size_size(hnd.begin(), ijk.begin()));
            thrust::transform(
              thrust::make_zip_iterator(make_tuple(y.begin(), j.begin(), C_fre,        C_hnd       )), // input - begin
              thrust::make_zip_iterator(make_tuple(y.end(),   j.end(),   C_fre+n_part, C_hnd+n_part)), // input - end
              y.begin(), // output
              detail::adve_helper<real_t>(opts_init.dy, apply)
            );
          }

          if (n_dims > 1) 
          {
            const pi_real_size
              C_abv(courant_z.begin(), pi_size_size(abv.begin(), ijk.begin())),
              C_blw(courant_z.begin(), pi_size_size(blw.begin(), ijk.begin()));
            thrust::transform(
              thrust::make_zip_iterator(make_tuple(z.begin(), k.begin(), C_blw,        C_abv       )), // input - begin
              thrust::make_zip_iterator(make_tuple(z.end(),   k.end(),   C_blw+n_part, C_abv+n_part)), // input - end
              z.begin(), // output
              detail::adve_helper<real_t>(opts_init.dz, apply)
            );
          }

          break; 
        }
        case 0: break;
        default: assert(false);
      }
    }

/*
    // apply change of position due to advection
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::adve_shift_apply()
    {   
      namespace arg = thrust::placeholders;
      switch (n_dims)
      {
        case 3:
        case 2:
        case 1:
        {
          thrust::transform( x.begin(), x.end(), dx.begin(), x.begin(), arg::_1 + arg::_2);
          if (n_dims > 2)
            thrust::transform( y.begin(), y.end(), dy.begin(), y.begin(), arg::_1 + arg::_2);
          if (n_dims > 1) 
            thrust::transform( z.begin(), z.end(), dz.begin(), z.begin(), arg::_1 + arg::_2);
          break; 
        }
        case 0: break;
        default: assert(false);
      }
    }
*/

    // calculate new positions using the predictor-corrector method with nearest-neighbour interpolation
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::adve()
    {   
      if(n_dims==0) return;
      namespace arg = thrust::placeholders;

      // old positions storage
      thrust_device::vector<real_t> 
        &x_old(tmp_device_real_part),
        &y_old(tmp_device_real_part2),
        &z_old(tmp_device_real_part1);

      // shift to coordiante system starting at halo's left edge
      thrust::transform(x.begin(), x.end(), x.begin(), arg::_1 + opts_init.dx);

      hskpng_ijk();        // get cell indices in new coordinates

      // save old x, y and z
      thrust::copy(x.begin(), x.end(), x_old.begin());
      if (n_dims > 2)
        thrust::copy(y.begin(), y.end(), y_old.begin());
      if (n_dims > 1) 
        thrust::copy(z.begin(), z.end(), z_old.begin());

      // ---- predictor step ----
      adve_calc(true);   // change x,y,z to mid-step values
      // due to numerics we could end up out of domain in z direction - move them back into domain since it would break next ijk
      if (n_dims > 1)
      {
        thrust::replace_if(
          z.begin(), z.end(), 
          arg::_1 >= opts_init.z1,     // condition
          opts_init.z1 - 1e-8 * opts_init.dz // TODO: sth smarter
        );
        thrust::replace_if(
          z.begin(), z.end(), 
          arg::_1 <= opts_init.z0,     // condition
          opts_init.z0 + 1e-8 * opts_init.dz // TODO: sth smarter
        );
      }

      // apply periodic boundary condition in y 
      if (n_dims == 3)
      {
        // adjust y_old to preserve y_old_post + y_1/2_bcnd = y_old_pre + y_1/2
        // TODO: do it in one call
        thrust::transform_if(
          y_old.begin(), y_old.end(), // input
          y.begin(), // strencil
          y_old.begin(), //out
          arg::_1 + (opts_init.y1 - opts_init.y0), // operation
          arg::_1 >= opts_init.y1 // condition
        );
        thrust::transform_if(
          y_old.begin(), y_old.end(), // input
          y.begin(), // strencil
          y_old.begin(), //out
          arg::_1 - (opts_init.y1 - opts_init.y0), // operation
          arg::_1 < opts_init.y0 // condition
        );
        //change y, copied from bcnd...
        thrust::transform(
          y.begin(), y.end(),
          y.begin(),
          detail::periodic<real_t>(opts_init.y0, opts_init.y1)
        );
      }
      // cell indices after predictor step
      hskpng_ijk();

      // ---- corrector step ----
      // save (x(t+1/2) + x(t)) in old x,y,z
      thrust::transform(
        x.begin(), x.end(), // arg1
        x_old.begin(),       // arg2
        x_old.begin(),       // out
        (arg::_1 + arg::_2)  // oper
      );
      if (n_dims > 2)
        thrust::transform(
          y.begin(), y.end(), // arg1
          y_old.begin(),       // arg2
          y_old.begin(),       // out
          (arg::_1 + arg::_2)  // oper
        );
      if (n_dims > 1) 
        thrust::transform(
          z.begin(), z.end(), // arg1
          z_old.begin(),       // arg2
          z_old.begin(),       // out
          (arg::_1 + arg::_2) // oper
        );

      // calculate rhs at midpoint position
      adve_calc(false);

      // calculate final position x(t+1) = (x(t+1/2) + x(t)) / 2 + 1/2 dx(x(t+1/2))
      thrust::transform(
        x.begin(), x.end(),  // arg1 - dx(x(t+1/2))
        x_old.begin(),       // arg2 - (x(t+1/2) + x(t))
        x.begin(),           // out
        (arg::_1 + arg::_2) / real_t(2.)       // oper
      );
      if (n_dims > 2)
        thrust::transform(
          y.begin(), y.end(), // arg1
          y_old.begin(),       // arg2
          y.begin(),       // out
          (arg::_1 + arg::_2) / real_t(2.)       // oper
        );
      if (n_dims > 1) 
        thrust::transform(
          z.begin(), z.end(), // arg1
          z_old.begin(),       // arg2
          z.begin(),       // out
          (arg::_1 + arg::_2) / real_t(2.)       // oper
        );
      // shift back to regular coordiante system
      thrust::transform(x.begin(), x.end(), x.begin(), arg::_1 - opts_init.dx);
    }
  };  
};
