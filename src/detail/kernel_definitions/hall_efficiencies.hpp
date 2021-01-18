#pragma once
namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class real_t>
      real_t hall_r_max()
      {
        return 300.e-6; // not used?
      }

      const int hall_n_rad = 11;
      const int hall_n_rat = 21;

      template<class real_t>
      void hall_efficiencies(std::vector<real_t> &vec)
      {
        //NVCC doesn't support extended initializer list. First make an array, then copy it to a vector.
        const real_t arr[hall_n_rat][hall_n_rad] = 
        {
          {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
          {.0001, .0001, .0001, .001, .005, .05, .2, .5, .77, .87, .97},
          {.0001, .0001, .002, .07, .4, .43, .58, .79, .93, .96, 1.},
          {.0001, .005, .02, .28, .6, .64, .75, .91, .97, .98, 1.},
          {.014, .016, .04, .5, .7, .77, .84, .95, .97, 1., 1.},
          {.017, .022, .085, .62, .78, .84, .88, .95, 1., 1., 1.},
          {.019, .03, .17, .68, .83, .87, .9, 1., 1., 1, 1},
          {.022, .043, .27, .74, .86, .89, .92, 1, 1, 1, 1},
          {.027, .05, .4, .78, .88, .9, .94, 1, 1, 1, 1},
          {.03, .064, .5, .8, .9, .91, .95, 1, 1, 1, 1},
          {.033, .072, .55, .8, .9, .91, .95, 1, 1, 1, 1},
          {.035, .079, .58, .8, .9, .91, .95, 1, 1, 1, 1},
          {.037, .082, .59, .78, .9, .91, .95, 1, 1, 1, 1},
          {.038, .08, .58, .77, .89, .91, .95, 1, 1, 1, 1},
          {.038, .076, .54, .76, .88, .92, .95, 1, 1, 1, 1},
          {.037, .067, .51, .77, .88, .93, .97, 1, 1, 1, 1},
          {.036, .057, .49, .77, .89, .95, 1, 1, 1, 1, 1},
          {.035, .048, .47, .78, .92, 1, 1.02, 1, 1, 1, 1},
          {.032, .04, .45, .79, 1.01, 1.03, 1.04, 1., 1., 1. ,1.},
          {.029, .033, .47, .95, 1.3, 1.7, 2.3, 1., 1., 1, 1},
          {.027, .027, .52, 1.4, 2.3, 3, 4, 1, 1, 1, 1}
        };

        vec = std::vector<real_t>(&(arr[0][0]), &(arr[hall_n_rat][hall_n_rad]));
      }

      template<class real_t>
      void hall_radii(std::vector<real_t> &vec)
      {
        //NVCC doesn't support extended initializer list. First make an array, then copy it to a vector.
        const real_t arr[hall_n_rad] = {10.e-6,20.e-6,30.e-6,40.e-6,50.e-6,60.e-6,70.e-6,100.e-6,150.e-6,200.e-6,300.e-6}; 
  
        vec = std::vector<real_t>(arr, arr + hall_n_rad );
      }
  
      template<class real_t>
      void hall_ratios(std::vector<real_t> &vec)
      {
        //NVCC doesn't support extended initializer list. First make an array, then copy it to a vector.
        const real_t arr[hall_n_rat] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5, 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
  
        vec = std::vector<real_t>(arr, arr + hall_n_rat );
      }
    }
  }
}
