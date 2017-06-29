#pragma once
namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class real_t>
      real_t hall_davis_no_waals_r_max()
      {
        return 300.e-6; // not used??
      }

      const int hall_davis_no_waals_n_rad = 15;
      const int hall_davis_no_waals_n_rat = 21;

      template<class real_t>
      void hall_davis_no_waals_efficiencies(std::vector<real_t> &vec)
      {
        //NVCC doesn't support extended initializer list. First make an array, then copy it to a vector.
        const real_t arr[hall_davis_no_waals_n_rat][hall_davis_no_waals_n_rad] = 
        {
          {0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0,
          0, 0, 0.001, 0.084, 0.05,
          0.2543, 0.5, 0.7807, 0.87, 0.97},
          {0, 0, 0, 0, 0,
          0, 0, 0.07, 0.4, 0.43,
          0.58, 0.79, 0.93, 0.96, 1},
          {0, 0, 0, 0, 0,
          0, 0, 0.28, 0.62, 0.64,
          0.7629, 0.91, 0.9687, 0.98, 1},
          {0, 0.01215, 0, 0, 0,
          0, 0, 0.5, 0.7, 0.77,
          0.84, 0.95, 0.95, 1, 1},
          {0.01122, 0.01215, 0.01742, 0.0209, 0.02843,
          0, 0.1357, 0.62, 0.79, 0.84,
          0.8829, 0.95, 1, 1, 1},
          {0.01122, 0.01215, 0.01742, 0.0261, 0.0341,
          0.08008, 0.1807, 0.68, 0.83, 0.87,
          0.9, 1, 1, 1, 1},
          {0.01122, 0.01641, 0.02219, 0.0261, 0.03975,
          0.09509, 0.2707, 0.74, 0.864, 0.89,
          0.9229, 1, 1, 1, 1},
          {0.01122, 0.01641, 0.02219, 0.03131, 0.04541,
          0.1101, 0.3157, 0.78, 0.88, 0.9,
          0.94, 1, 1, 1, 1},
          {0.01503, 0.02067, 0.02696, 0.03651, 0.05107,
          0.1251, 0.4057, 0.8, 0.9, 0.91,
          0.95, 1, 1, 1, 1},
          {0.01503, 0.02067, 0.02696, 0.03822, 0.05673,
          0.143, 0.4507, 0.8, 0.9, 0.91,
          0.95, 1, 1, 1, 1},
          {0.01503, 0.02067, 0.02696, 0.03822, 0.05381,
          0.1336, 0.39, 0.8, 0.9, 0.91,
          0.95, 1, 1, 1, 1},
          {0.02, 0.02264, 0.02696, 0.03644, 0.05089,
          0.1242, 0.3597, 0.78, 0.9, 0.91,
          0.95, 1, 1, 1, 1},
          {0.02, 0.02264, 0.02696, 0.03466, 0.04797,
          0.1148, 0.299, 0.77, 0.888, 0.91,
          0.95, 1, 1, 1, 1},
          {0.02, 0.0246, 0.02696, 0.03288, 0.04504,
          0.09603, 0.2687, 0.76, 0.88, 0.92,
          0.95, 1, 1, 1, 1},
          {0.02497, 0.0246, 0.02696, 0.03288, 0.04212,
          0.08665, 0.2081, 0.77, 0.882, 0.93,
          0.9743, 1, 1, 1, 1},
          {0.02497, 0.0246, 0.02696, 0.0311, 0.0392,
          0.07726, 0.1778, 0.77, 0.89, 0.95,
          1, 1, 1, 1, 1},
          {0.02497, 0.02657, 0.02696, 0.02932, 0.03627,
          0.06788, 0.1171, 0.78, 0.938, 1,
          1.023, 1, 1, 1, 1},
          {0.02497, 0.02657, 0.02696, 0.02843, 0.03335,
          0.0538, 0.08679, 0.79, 1.01, 1.03,
          1.04, 1, 1, 1, 1},
          {0.02696, 0.02696, 0.02696, 0.02843, 0.03335,
          0.0538, 0.08679, 0.95, 1.5, 1.7,
          2.543, 1, 1, 1, 1},
          {0.02696, 0.02696, 0.02696, 0.02843, 0.03335,
          0.0538, 0.08679, 1.4, 2.3, 3,
          4, 1, 1, 1, 1}
        };

        vec = std::vector<real_t>(&(arr[0][0]), &(arr[hall_davis_no_waals_n_rat][hall_davis_no_waals_n_rad]));
      }

      template<class real_t>
      void hall_davis_no_waals_radii(std::vector<real_t> &vec)
      {
        //NVCC doesn't support extended initializer list. First make an array, then copy it to a vector.
        const real_t arr[hall_davis_no_waals_n_rad] = {6.e-6,8.e-6,10.e-6,15.e-6,20.e-6,25.e-6,30.e-6,40.e-6,50.e-6,60.e-6,70.e-6,100.e-6,150.e-6,200.e-6,300.e-6}; 
  
        vec = std::vector<real_t>(arr, arr + hall_davis_no_waals_n_rad );
      }
  
      template<class real_t>
      void hall_davis_no_waals_ratios(std::vector<real_t> &vec)
      {
        //NVCC doesn't support extended initializer list. First make an array, then copy it to a vector.
        const real_t arr[hall_davis_no_waals_n_rat] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5, 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
  
        vec = std::vector<real_t>(arr, arr + hall_davis_no_waals_n_rat );
      }
    }
  }
}
