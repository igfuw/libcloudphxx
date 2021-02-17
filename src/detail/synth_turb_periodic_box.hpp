#pragma once
#include "urand.hpp"
#include "../../include/libcloudph++/lgrngn/backend.hpp"

namespace libcloudphxx
{
  namespace synth_turb_periodic_box
  {
    // TODO: a co z vectors.size==0?
    inline int mode_degeneracy_generator(const int &n, std::vector<std::array<int,3>> &vectors) // returns degeneracy, fills vectors with half of the degeneracy (opposite vectors are not added)
    {
      int degeneracy = 0;
      vectors.clear();
      const double max = sqrt(n);
      for(int nx=-max; nx<=max; ++nx)
        for(int ny=-max; ny<=max; ++ny)
          for(int nz=-max; nz<=max; ++nz)
          {
            if(nx*nx + ny*ny + nz*nz == n)
            {
   //           std::cout << "(" << nx << "," << ny << "," << nz << ")" << std::endl;
              ++degeneracy;
              if(std::find(vectors.begin(), vectors.end(), std::array<int,3>{-nx,-ny,-nz}) == vectors.end()) // add only if the opposite vector has not been added yet
                vectors.push_back({nx,ny,nz});
            }
          }
   //   std::cout << "degeneracy: " << degeneracy << std::endl;
      return degeneracy;
    }

    template <class real_t>
    class mode
    {
      private:
      real_t std_dev, // sqrt(variance)
             wn;      // frequency

      unsigned int Nwaves; // actual number of wave vectors in thath mode
      // wavevectors in the form k = (nx,ny,nz) * 2 PI / L, where n is integer to get periodic flow
//      const int nn; // = nx^2 + ny^2 + nz^2

      // arrays of size [3][Nwaves]
      thrust_device::vector<real_t> enm, knm, Anm, Bnm;

      public:
      BOOST_GPU_ENABLED
      void update_time(const real_t &dt)
      {
        /*
        std::normal_distribution<real_t> normal_d(0,1);
        std::default_random_engine local_rand_eng(std::random_device{}());
        real_t relax = exp(-wn * dt);

        for(int m=0; m<Nwaves; m+=2)
        {
          for(int i=0; i<3; ++i)
          {
            Anm[i][m] = relax * Anm[i][m] + std_dev * sqrt(1. - relax * relax) * normal_d(local_rand_eng);
            Anm[i][m+1] = -Anm[i][m];

            Bnm[i][m] = relax * Bnm[i][m] + std_dev * sqrt(1. - relax * relax) * normal_d(local_rand_eng);
            Bnm[i][m+1] = Bnm[i][m];
          }
        }
        */
      }

      BOOST_GPU_ENABLED
      mode(const real_t &std_dev,
           const real_t &wn,
           const int Nwaves,
           const std::vector<real_t> enm,
           const std::vector<real_t> Anm,
           const std::vector<real_t> Bnm,
           const std::vector<real_t> knm
           ):
        Nwaves(Nwaves),
        std_dev(std_dev),
        wn(wn),
        enm(enm),
        Anm(Anm),
        Bnm(Bnm),
        knm(knm)
      {
      }
    };

    template <class real_t, lgrngn::backend_t device, int Nmodes, int Nwaves_max>
    class synth_turb
    {
      private:
      thrust_device::vector<mode<real_t>> modes;
      thrust_device::vector<real_t> &tmp_device_real_part;
      lgrngn::detail::rng<real_t, device> rng; // separate rng, because we want the same sequence of random numbers on all distmem GPUS/nodes

      void update_time(const real_t &dt)
      {
        // get random normally distributed vars
        int total_number_of_modes = 0;
        rng.generate_normal_n(tmp_device_real_part, total_number_of_modes);

        for(int n=0; n<Nmodes; ++n)
        {
        }
      }

      public:

      synth_turb(
        const real_t &eps,        // TKE dissipation rate [m2/s3], for now it has to be a constant value
        int seed,                 // rng seed
        thrust_device::vector<real_t> &_tmp_device_real_part, // temp space for real_t on the device
        const real_t &Lmax = 100, // maximum length scale [m]
        const real_t &Lmin = 1e-3 // Kolmogorov length scale [m]
      ):
        rng(seed),
        tmp_device_real_part(_tmp_device_real_part)
      {
        if(Nwaves_max % 2 != 0) throw std::runtime_error("Nwaves_max needs to be even, because we need to include opposites of all wavevectors.");

        const real_t lambda = 1; // unsteadiness parameter, within [0,1]; see Sidin et al. 2009

        real_t k[Nmodes],   // norms of wave vectors
               dk[Nmodes],  // differences between norms of subsequent wave vectors
               E[Nmodes],   // kinetic energy in a mode
               std_dev[Nmodes], // sqrt(variances)
               wn[Nmodes];   // frequencies

        int nn[Nmodes],     // nn = nx^2 + ny^2 + nz^2
            Nwaves[Nmodes]; // actual number of wave vectors in thath mode

        std::array<std::vector<real_t>, Nmodes> enm, // unit vectors along wavevectors
                                                Anm, Bnm, // coefficients of the Fourier transform
                                                knm;

  
        // --- linear distribution of nn (nn = 1, 2, 3, 4, ..., Nmodes) ---
        /*
        {
          for(int n=0; n<Nmodes; ++n)
            nn[n] = n+1;
        }
        */
        // --- geometric distribution of nn ---
        {
          if(Nmodes > Lmax / Lmin)
            throw std::runtime_error("too many modes: Nmodes is greater than Lmax / Lmin");
  
          nn[0]=1;
  
          real_t alpha = pow(Lmax / Lmin, 1. / (Nmodes - 1));
          while(1)
          {
            for(int n=1; n<Nmodes; ++n)
            {
              nn[n] = -1;
              int exponent = n;
              while(nn[n] <= nn[n-1])
              {
                nn[n] = std::round(std::pow(alpha, exponent++));
              }
    //          std::cerr << "alpha: " << alpha << " nn[" << n << "]: " << nn[n] << std::endl;
              if(nn[n] > Lmax / Lmin) break;
            }
            if(nn[Nmodes-1] <= Lmax / Lmin && nn[Nmodes-1]!=0)
              break;
            else
              alpha /= 1.001;
          }
          assert(nn[Nmodes-1] > 0);
        }

        for(int n=0; n<Nmodes; ++n)
        {
  //        std::cerr << "nn[" << n << "]: " << nn[n] << std::endl;
          k[n] = sqrt(real_t(nn[n])) * (2. * M_PI / Lmax);
        //  k[n] = sqrt(nn[n]) * (2. * M_PI / Lmax);
        }
  
        std::vector<std::array<int,3>> vectors;
        for(int n=0; n<Nmodes; ++n)
        {
          Nwaves[n] = mode_degeneracy_generator(nn[n], vectors);
  
          if(Nwaves[n] > Nwaves_max) // random shuffle, because not all possible degeneracies will be used
          {
            std::default_random_engine local_rand_eng(std::random_device{}());
            std::shuffle(std::begin(vectors), std::end(vectors), local_rand_eng);
            Nwaves[n] = Nwaves_max;
          }
  //        if(Nwaves_max != 6) throw std::runtime_error("nwaves max needs to be 6 for this test");
  //        vectors = {{1,0,0},{0,1,0},{0,0,1}};
  
          for(int m=0; m<Nwaves[n]; m+=2)
          {
          //  enm[n][0][m] = vectors.at(m/2)[0] / sqrt(real_t(nn[n]));
          //  enm[n][1][m] = vectors.at(m/2)[1] / sqrt(real_t(nn[n]));
          //  enm[n][2][m] = vectors.at(m/2)[2] / sqrt(real_t(nn[n]));
          //  // opposite vector
          //  enm[n][0][m+1] = -vectors.at(m/2)[0] / sqrt(real_t(nn[n]));
          //  enm[n][1][m+1] = -vectors.at(m/2)[1] / sqrt(real_t(nn[n]));
          //  enm[n][2][m+1] = -vectors.at(m/2)[2] / sqrt(real_t(nn[n]));

            enm[n].push_back(vectors.at(m/2)[0] / sqrt(real_t(nn[n])));
            enm[n].push_back(vectors.at(m/2)[1] / sqrt(real_t(nn[n])));
            enm[n].push_back(vectors.at(m/2)[2] / sqrt(real_t(nn[n])));
            // opposite vector
            enm[n].push_back(-vectors.at(m/2)[0] / sqrt(real_t(nn[n])));
            enm[n].push_back(-vectors.at(m/2)[1] / sqrt(real_t(nn[n])));
            enm[n].push_back(-vectors.at(m/2)[2] / sqrt(real_t(nn[n])));
          }
        }
  
        // Energy spectrum; first equation in Appendix of Sidin et al. 2009, but without the corrections f_L and f_eta
        for(int n=0; n<Nmodes; ++n)
          E[n] = 1.44 * pow(eps, 2./3.) * pow(k[n], -5./3.);
  
        // Wave vector differences
        dk[0] =         0.5 *  (k[1]        - k[0]);
        dk[Nmodes-1] =  0.5 *  (k[Nmodes-1] - k[Nmodes-2]);
  
        for(int n=1; n<Nmodes-1; ++n)
          dk[n] = 0.5 * (k[n+1] - k[n-1]);

        // std deviations
        for(int n=0; n<Nmodes; ++n)
        {
          std_dev[n] = sqrt(E[n] * dk[n] / Nwaves[n]);
        }
  
        // Frequencies; Eq. A7 in Sidin et al. 2009
        for(int n=0; n<Nmodes; ++n)
          wn[n] = lambda * sqrt(E[n] * k[n] * k[n] * k[n]);
  
        // std::cerr << "SynthTurb3d ctor debug output:" << std::endl;
        // for(int n=0; n<Nmodes; ++n)
        // {
        //    std::cerr << "n: " << n << " k[n]: " << k[n] << std::endl;
        //    std::cerr << "n: " << n << " E[n]: " << E[n] << std::endl;
        //    std::cerr << "n: " << n << " dk[n]: " << dk[n] << std::endl;
        //    std::cerr << "n: " << n << " wn[n]: " << wn[n] << std::endl;
        //    std::cerr << std::endl;
        // }
  
        std::default_random_engine local_rand_eng(std::random_device{}());
        for(int n=0; n<Nmodes; ++n)
        {
          std::normal_distribution<real_t> G_d(0, std_dev[n]);
  
          for(int m=0; m<Nwaves[n]; ++m)
          {
            // knm = unit vector * magnitude
            for(int i=0; i<3; ++i)
           //   knm[n][i][m] = enm[n][i][m] * k[n];
              knm[n].push_back(enm[n].at(3*m+i) * k[n]);
  
            // init random coefficients

            /*
            for(int i=0; i<3; ++i)
            {
              Anm[i][n][m] = 0;
              Bnm[i][n][m] = 0;
            }
  */
  
            for(int i=0; i<3; ++i)
            {
           //   Anm[n][i][m] = G_d(local_rand_eng);
           //   Bnm[n][i][m] = G_d(local_rand_eng);
              Anm[n].push_back(G_d(local_rand_eng));
              Bnm[n].push_back(G_d(local_rand_eng));
            }
          }
        }

        // init modes
        for(int n=0; n<Nmodes; ++n)
        {
          modes.push_back(mode<real_t>(
            std_dev[n],
            wn[n],
            Nwaves[n],
            enm[n],
            Anm[n],
            Bnm[n],
            knm[n]
          ));
          /*
      mode(const real_t &std_dev,
           const real_t &wn,
           const int &Nwaves,
           const real_t enm[3][Nwaves_max],
           const real_t Anm[3][Nwaves_max],
           const real_t Bnm[3][Nwaves_max],
           const real_t knm[3][Nwaves_max]
           )
      */
        }
      }
    };

    // TODO: 3D, 2D, 1D
  };
};


                            
