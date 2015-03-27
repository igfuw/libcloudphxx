#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, typename n_t>
    BOOST_GPU_ENABLED
    real_t kernel_geometric<real_t, n_t>::bilinear_interpolation(real_t r1, real_t r2) const
    {
      n_t dx, dy, // distance between efficiencies in the matrix
          x[4];   // positions in the (R,r) space of the defined efficiencies. x1, x2, y1, y2

      if(r1 >= kernel_base<real_t, n_t>::r_max) 
        r1 = kernel_base<real_t, n_t>::r_max - 1e-6; //TODO : improve this
      if(r2 >= kernel_base<real_t, n_t>::r_max) 
        r2 = kernel_base<real_t, n_t>::r_max - 1e-6; //TODO : improve this

      if(r1 >= 100.)
      {
        x[0] = int(r1/10.) * 10; //TODO : what if r1 exactly equal to, for example, 120.?
        dx = 10;
      }
      else
      {
        x[0] = int(r1);
        dx = 1;
      }

      if(r2 >= 100.)
      {
        x[2] = int(r2/10.) * 10;
        dy = 10;
      }
      else
      {
        x[2] = int(r2);
        dy = 1;
      }
      x[1] = x[0] + dx;
      x[3] = x[2] + dy;

      n_t iv[4];     // kernel_parameters vector indices of the four neighbouring efficiencies

      //kernel_index does the same type casting as just done above?
      iv[0] = detail::kernel_vector_index<n_t>(detail::kernel_index<real_t, n_t>(x[0]), detail::kernel_index<real_t, n_t>(x[2]));
      iv[1] = detail::kernel_vector_index<n_t>(detail::kernel_index<real_t, n_t>(x[1]), detail::kernel_index<real_t, n_t>(x[2]));
      iv[2] = detail::kernel_vector_index<n_t>(detail::kernel_index<real_t, n_t>(x[0]), detail::kernel_index<real_t, n_t>(x[3]));
      iv[3] = detail::kernel_vector_index<n_t>(detail::kernel_index<real_t, n_t>(x[1]), detail::kernel_index<real_t, n_t>(x[3]));

      real_t w[4];   //  weighting factors
      w[0] = r1 - x[0];
      w[1] = x[1] - r1;
      w[2] = r2 - x[2];
      w[3] = x[3] - r2;

      return
      (
        kernel_base<real_t, n_t>::k_params[iv[0]] * w[1] * w[3] +
        kernel_base<real_t, n_t>::k_params[iv[1]] * w[0] * w[3] +
        kernel_base<real_t, n_t>::k_params[iv[2]] * w[1] * w[2] +
        kernel_base<real_t, n_t>::k_params[iv[3]] * w[0] * w[2]
      ) / dx / dy;
    }
  }
}
