#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    // a version for Jon's efficiency tables
    template <typename real_t, typename n_t>
    BOOST_GPU_ENABLED
    real_t kernel_with_efficiencies<real_t, n_t>::interpolated_efficiency(real_t r1, real_t r2) const //radii in meters
    {
      r1*=1e6; r2*=1e6; // to work on micrometers

      if(r1 >= r_max) 
        r1 = r_max - 1e-6;
      if(r2 >= r_max) 
        r2 = r_max - 1e-6; 

      n_t dx, dy, // distance between efficiencies in the matrix
          x[4];   // positions in the (R,r) space of the defined efficiencies. x1, x2, y1, y2

      if(r1 >= 100.)
      {
        x[0] = floor(r1/10.) * 10;
        dx = 10;
      }
      else
      {
        x[0] = floor(r1);
        dx = 1;
      }

      if(r2 >= 100.)
      {
        x[2] = floor(r2/10.) * 10;
        dy = 10;
      }
      else
      {
        x[2] = floor(r2);
        dy = 1;
      }
      x[1] = x[0] + dx;
      x[3] = x[2] + dy;

      thrust_size_t iv[4];     // kernel_parameters vector indices of the four neighbouring efficiencies

      iv[0] = detail::kernel_vector_index<n_t>(detail::kernel_index<n_t>(x[0]), detail::kernel_index<n_t>(x[2]));
      iv[1] = detail::kernel_vector_index<n_t>(detail::kernel_index<n_t>(x[1]), detail::kernel_index<n_t>(x[2]));
      iv[2] = detail::kernel_vector_index<n_t>(detail::kernel_index<n_t>(x[0]), detail::kernel_index<n_t>(x[3]));
      iv[3] = detail::kernel_vector_index<n_t>(detail::kernel_index<n_t>(x[1]), detail::kernel_index<n_t>(x[3]));

      real_t w[4];   //  weighting factors
      w[0] = r1 - x[0];
      w[1] = x[1] - r1;
      w[2] = r2 - x[2];
      w[3] = x[3] - r2;

      return
      (
        k_coll_eff[iv[0]] * w[1] * w[3] +
        k_coll_eff[iv[1]] * w[0] * w[3] +
        k_coll_eff[iv[2]] * w[1] * w[2] +
        k_coll_eff[iv[3]] * w[0] * w[2]
      ) / dx / dy;
    }

    // a version for efficiency tables with defined radius and ratio values
    template <typename real_t, typename n_t>
    BOOST_GPU_ENABLED
    real_t kernel_with_efficiencies<real_t, n_t>::interpolated_efficiency_radrat(real_t r1, real_t r2) const //radii in meters
    {
      if(r1==r2) return 0.;
#if !defined(__NVCC__)
      using std::max;
      using std::min;
#endif
      const real_t radius = max(r1,r2);
      const real_t ratio = min(r1,r2) / radius;
      real_t efficiency;

      // find indexes of first not less value of radius and ratio
      int rad_pos = n_rad, rat_pos; 
      for(int i = 0; i < n_rad; ++i) if(k_coll_eff_rad[i] > radius) {rad_pos = i; break;} // within [0, n_rad]
      for(int i = 1; i < n_rat; ++i) if(k_coll_eff_rat[i] > ratio)  {rat_pos = i; break;} // within [1, n_rat], since k_coll_eff_rat should encompass the whole [0,1] range

      if(rad_pos < n_rad)
      {
        if(rad_pos >= 1)
        {
          real_t p = (radius - k_coll_eff_rad[rad_pos-1]) / (k_coll_eff_rad[rad_pos] - k_coll_eff_rad[rad_pos - 1]);
          real_t q = (ratio - k_coll_eff_rat[rat_pos-1]) / (k_coll_eff_rat[rat_pos] - k_coll_eff_rat[rat_pos - 1]);
          efficiency = (1-p)*(1-q)*k_coll_eff[(rat_pos-1)*n_rad + rad_pos-1] +
                       p*(1-q)*    k_coll_eff[(rat_pos-1)*n_rad + rad_pos  ] +
                       q*(1-p)*    k_coll_eff[rat_pos*n_rad     + rad_pos-1] +
                       p*q*        k_coll_eff[rat_pos*n_rad     + rad_pos  ];
        }
        else
        {
          real_t q = (ratio - k_coll_eff_rat[rat_pos-1]) / (k_coll_eff_rat[rat_pos] - k_coll_eff_rat[rat_pos - 1]);
          efficiency = (1-q)*k_coll_eff[(rat_pos-1) * n_rad] + q*k_coll_eff[rat_pos*n_rad];
        }
      }
      else
      {
        real_t q = (ratio - k_coll_eff_rat[rat_pos-1]) / (k_coll_eff_rat[rat_pos] - k_coll_eff_rat[rat_pos - 1]);
        efficiency = (1-q)*k_coll_eff[rat_pos*n_rad-1] + q*k_coll_eff[(rat_pos+1)*n_rad-1];
        efficiency = efficiency > 1.? 1. : efficiency;
      }
      return efficiency;
    }
  }
}
