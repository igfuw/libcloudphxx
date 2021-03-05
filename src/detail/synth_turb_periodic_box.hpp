#pragma once
#include <random>
#include "urand.hpp"
#include "thrust.hpp"
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

      // pointers to arrays of size [3][Nwaves] stored in device memory
      thrust_device::pointer<real_t> enm, knm, Anm, Bnm;

      public:
      BOOST_GPU_ENABLED
      void update_time(const real_t &dt, const real_t &rand_normal)
      {
        const real_t relax = exp(-wn * dt);
        const real_t *p_rand_normal(&rand_normal); // we rely on iterators passing pointers do data, without copies!

        for(int m=0; m<Nwaves; m+=2)
        {
          for(int i=0; i<3; ++i)
          {
            Anm[3*m+i] = relax * Anm[3*m+i] + std_dev * sqrt(1. - relax * relax) * p_rand_normal[3*m+2*i];
            Anm[3*(m+1)+i] = -Anm[3*m+i];

            Bnm[3*m+i] = relax * Bnm[3*m+i] + std_dev * sqrt(1. - relax * relax) * p_rand_normal[3*m+2*i+1];
            Bnm[3*(m+1)+i] = Bnm[3*m+i];
          }
        }
      }

      BOOST_GPU_ENABLED
      void add_vel(const real_t &x, const real_t &y, const real_t &z, real_t &u, real_t &v, real_t &w)
      {
        for(int m=0; m<Nwaves; ++m)
        {
          const real_t r = knm[3*m] * x + knm[3*m+1] * y + knm[3*m+2] * z;
          const real_t cosr = cos(r);
          const real_t sinr = sin(r);

          u += (Anm[3*m+1] * enm[3*m+2] - Anm[3*m+2] * enm[3*m+1])*cosr - (Bnm[3*m+1] * enm[3*m+2] - Bnm[3*m+2] * enm[3*m+1])*sinr;
          v += (Anm[3*m+2] * enm[3*m]   - Anm[3*m]   * enm[3*m+2])*cosr - (Bnm[3*m+2] * enm[3*m]   - Bnm[3*m]   * enm[3*m+2])*sinr;
          w += (Anm[3*m]   * enm[3*m+1] - Anm[3*m+1] * enm[3*m])  *cosr - (Bnm[3*m]   * enm[3*m+1] - Bnm[3*m+1] * enm[3*m])*sinr;
        }
      }

      BOOST_GPU_ENABLED
      mode(const real_t &std_dev,
           const real_t &wn,
           const int Nwaves,
           const thrust_device::pointer<real_t> enm,
           const thrust_device::pointer<real_t> Anm,
           const thrust_device::pointer<real_t> Bnm,
           const thrust_device::pointer<real_t> knm
           ):
        Nwaves(Nwaves),
        std_dev(std_dev),
        wn(wn),
        enm(enm),
        Anm(Anm),
        Bnm(Bnm),
        knm(knm)
      {}
    };

    template<class real_t>
    struct mode_update_time
    {
      const real_t dt;

      mode_update_time(const real_t &dt): dt(dt) {}

      template <class tpl_t>
      BOOST_GPU_ENABLED
      void operator()(const tpl_t &tpl)
      {
        thrust::get<0>(tpl).update_time(dt, thrust::get<1>(tpl));
      }
    };

    template<class real_t>
    struct mode_add_vel
    {
      thrust_device::pointer<mode<real_t>> modep;

      mode_add_vel(thrust_device::pointer<mode<real_t>> modep): modep(modep) {}

      template <class tpl_t>
      BOOST_GPU_ENABLED
      void operator()(const tpl_t &tpl_pos_vel) //const mode<real_t> *mp, thrust_device::pointer<real_t> rand_normal)
      {
        thrust::raw_pointer_cast(modep)->add_vel(thrust::get<0>(tpl_pos_vel), thrust::get<1>(tpl_pos_vel), thrust::get<2>(tpl_pos_vel), thrust::get<3>(tpl_pos_vel), thrust::get<4>(tpl_pos_vel), thrust::get<5>(tpl_pos_vel));
      }
    };

    template <class real_t, lgrngn::backend_t device>
    class synth_turb
    {
      private:
      const int Nmodes, Nwaves_max;
      thrust_device::vector<mode<real_t>> modes;
      thrust_device::vector<real_t> &tmp_device_real_part;
      lgrngn::detail::rng<real_t, device> rng; // separate rng, because we want the same sequence of random numbers on all distmem GPUS/nodes

      std::vector<thrust_device::vector<real_t>> enm, // unit vectors along wavevectors
                                                 Anm, Bnm, // coefficients of the Fourier transform
                                                 knm;

      public:

      void update_time(const real_t &dt)
      {
        // get random normally distributed vars
        auto &rand_nrml(tmp_device_real_part);
        auto total_number_of_modes = modes.size();
        assert(total_number_of_modes * Nwaves_max * 6 <= rand_nrml.size());
        rng.generate_normal_n(rand_nrml, total_number_of_modes * Nwaves_max * 6); // each mode needs nwaves*6 random numbers; TODO: use actual nwaves not Nwaves_max

        namespace arg = thrust::placeholders;

        auto zip_mode_rand = 
          thrust::make_zip_iterator(thrust::make_tuple(
            modes.begin(),
            thrust::make_permutation_iterator(
              rand_nrml.begin(),
              thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                arg::_1 * Nwaves_max * 6
              )
            )
          ));

        thrust::for_each(
          zip_mode_rand,
          zip_mode_rand + modes.size(),
          mode_update_time<real_t>(dt)
        );

        for(int i=0; i<Nmodes; ++i)
        {
          std::cerr << "mode " << i << " Anm: " << std::endl;
          lgrngn::debug::print(Anm[i]);
          std::cerr << "mode " << i << " Bnm: " << std::endl;
          lgrngn::debug::print(Bnm[i]);
        }
      }

      template <class tpl_t, class thrust_size_t>
      void calc_vel(const tpl_t &zip_pos_vel, const thrust_size_t n_part)
      {
        for(int i=0; i<modes.size(); ++i)
        {
          thrust::for_each(
            zip_pos_vel,
            zip_pos_vel + n_part,
            mode_add_vel<real_t>(modes.data() + i)
          );
        }
      }

      synth_turb(
        const real_t &eps,        // TKE dissipation rate [m2/s3], for now it has to be a constant value
        const int &_Nmodes,
        const int &_Nwaves_max,
        const real_t &Lmax, // maximum length scale [m]
        const real_t &Lmin, // Kolmogorov length scale [m]
        const int &seed,                 // rng seed
        thrust_device::vector<real_t> &_tmp_device_real_part // temp space for real_t on the device
      ):
        Nmodes(_Nmodes),
        Nwaves_max(_Nwaves_max),
        rng(seed),
        tmp_device_real_part(_tmp_device_real_part)
      {
        if(Nwaves_max % 2 != 0) throw std::runtime_error("Nwaves_max needs to be even, because we need to include opposites of all wavevectors.");

        Anm.resize(Nmodes);
        Bnm.resize(Nmodes);
        knm.resize(Nmodes);
        enm.resize(Nmodes);

        const real_t lambda = 1; // unsteadiness parameter, within [0,1]; see Sidin et al. 2009

        std::vector<real_t> k(Nmodes),   // norms of wave vectors
                            dk(Nmodes),  // differences between norms of subsequent wave vectors
                            E(Nmodes),   // kinetic energy in a mode
                            std_dev(Nmodes), // sqrt(variances)
                            wn(Nmodes);   // frequencies

        std::vector<int> nn(Nmodes),     // nn = nx^2 + ny^2 + nz^2
                         Nwaves(Nmodes); // actual number of wave vectors in thath mode

  
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

std::cerr << "ST ctor mode " << n << std::endl;
std::cerr << "ST ctor nn " << nn[n] << std::endl;
std::cerr << "ST ctor Nwaves " << Nwaves[n] << std::endl;
  
          for(int m=0; m<Nwaves[n]; ++m)
          {
            // knm = unit vector * magnitude
            for(int i=0; i<3; ++i)
              knm[n].push_back(enm[n][3*m+i] * k[n]);
  
            // init random coefficients
  
            for(int i=0; i<3; ++i)
            {
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
            enm[n].data(),
            Anm[n].data(),
            Bnm[n].data(),
            knm[n].data()
          ));
        }
      }
    };

    // TODO: 3D, 2D, 1D
  };
};


                            
