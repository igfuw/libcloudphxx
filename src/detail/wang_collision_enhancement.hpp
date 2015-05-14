#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      // calculate collision efficiency enhancement due to turbulence
      // cf. Wang et al. (2009) atoms. sci. let. table 1
      template<class real_t>
      struct wang_collision_enhancement_t
      {
        static const real_t aR0[7];            // collector radius
        static const real_t arat[11];          // ratio of radii
        static const real_t aeta_e[11][2][7];  // [rat][eps][R0]

        //TODO: is it really more efficient to store it in device memory?
        thrust_device::vector<real_t> vR0;
        thrust_device::vector<real_t> vrat;
        thrust_device::vector<real_t> veta_e;

        thrust_device::pointer<real_t> R0;
        thrust_device::pointer<real_t> rat;
        thrust_device::pointer<real_t> eta_e;

        // ctor
        wang_collision_enhancement_t(): 
        vR0(aR0, aR0 + sizeof(aR0) / sizeof(aR0[0])), 
        vrat(arat, arat + sizeof(arat) / sizeof(arat[0])), 
        veta_e(&aeta_e[0][0][0], &aeta_e[0][0][0] + sizeof(aeta_e) / sizeof(aeta_e[0][0][0])),
        R0(vR0.data()),
        rat(vrat.data()),
        eta_e(veta_e.data())
        {}

        BOOST_GPU_ENABLED
        int eta_e_index(const int &n_rat, const int &n_eps, const int &n_R0) const  {return (7 * (2 * n_rat +  n_eps) + n_R0);}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &r1, const real_t &r2, const real_t &eps) const // radii in meters, dissipation rate in m^2/s^3
        {
          // find respective array positions
          real_t r, R;
          r1 > r2 ? (R = r1, r = r2) : (R = r2, r = r1);   // R > r
          if(R > 100e-6) return 1.;                    // no enhancement for large collector droplets
  
          int n_eps, n_R0, n_rat;
  
          eps <= 2.5e-2 ? n_eps = 0 : n_eps = 1;       // TODO: interpolate in epsilon space? cf. Benmoshe et al. 2012
  
          for(n_R0 = 0; n_R0 < 7; ++n_R0)
            if(R0[n_R0] > R) break;
        
          real_t ratio = r / R; 
          for(n_rat = 1; n_rat < 11; ++n_rat)
            if(rat[n_rat] > ratio) break;
  
          if(n_R0 == 0) return eta_e[eta_e_index(n_rat, n_eps, n_R0)]; // do not interpolate if both are smaller than 10um
  
          // two-dimensional linear interpolation
          // note: we interpolate in (radius, ratio) space, not in (radius, radius)...
          // note2: stagnant efficiency is also interpolated, hence we do it twice
          // TODO: prepare an already interpolated table of enhancement factors
          //       with exactly same layout as table of efficiencies, then
          //       first apply enhancement, interpolate afterward
  
          real_t w[4];   //  weighting factors
          w[0] = R - R0[n_R0-1];
          w[1] = R0[n_R0] - R; 
          w[2] = ratio - rat[n_rat-1];
          w[3] = rat[n_rat] - ratio;   
  
          return
          (   
            eta_e[eta_e_index(n_rat-1, n_eps, n_R0-1)] * w[1] * w[3] +
            eta_e[eta_e_index(n_rat-1, n_eps, n_R0)]   * w[0] * w[3] +
            eta_e[eta_e_index(n_rat, n_eps, n_R0-1)]   * w[1] * w[2] +
            eta_e[eta_e_index(n_rat, n_eps, n_R0)]     * w[0] * w[2]
          ) / (R0[n_R0] - R0[n_R0-1])
            / (rat[n_rat]-rat[n_rat-1]); 
        }
      };

     //definitions of member arrays
     template<class real_t>
     const real_t wang_collision_enhancement_t<real_t>::aR0[7] = {10e-6, 20e-6, 30e-6, 40e-6, 50e-6, 60e-6, 100e-6};
     template<class real_t>
     const real_t wang_collision_enhancement_t<real_t>::arat[11] = {0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.};
     template<class real_t>
     const real_t wang_collision_enhancement_t<real_t>::aeta_e[11][2][7] =  
     {                                       
       {{1.74, 1.74, 1.773,1.49, 1.207,1.207,1.0},
       {4.976,4.976,3.593,2.519,1.445,1.445,1.0},},
       {{1.46, 1.46, 1.421,1.245,1.069,1.069,1.0},
       {2.984,2.984,2.181,1.691,1.201,1.201,1.0},},
       {{1.32, 1.32, 1.245,1.123,1.000,1.000,1.0},
       {1.988,1.988,1.475,1.313,1.150,1.150,1.0},},
       {{1.250,1.250,1.148,1.087,1.025,1.025,1.0},
       {1.490,1.490,1.187,1.156,1.126,1.126,1.0},},
       {{1.186,1.186,1.066,1.060,1.056,1.056,1.0},
       {1.249,1.249,1.088,1.090,1.092,1.092,1.0},},
       {{1.045,1.045,1.000,1.014,1.028,1.028,1.0},
       {1.139,1.139,1.130,1.091,1.051,1.051,1.0},},
       {{1.070,1.070,1.030,1.038,1.046,1.046,1.0},
       {1.220,1.220,1.190,1.138,1.086,1.086,1.0},},
       {{1.000,1.000,1.054,1.042,1.029,1.029,1.0},
       {1.325,1.325,1.267,1.165,1.063,1.063,1.0},},
       {{1.223,1.223,1.117,1.069,1.021,1.021,1.0},
       {1.716,1.716,1.345,1.223,1.100,1.100,1.0},},
       {{1.570,1.570,1.244,1.166,1.088,1.088,1.0},
       {3.788,3.788,1.501,1.311,1.120,1.120,1.0},},
       {{20.3, 20.3, 14.6, 8.61, 2.60, 2.60, 1.0},
       {36.52,36.52,19.16,22.80,26.0, 26.0, 1.0}}
     };
    };
  };
};
