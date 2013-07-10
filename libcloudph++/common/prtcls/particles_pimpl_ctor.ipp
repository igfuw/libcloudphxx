// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      // pimpl stuff
      template <typename real_t, int device>
      struct particles<real_t, device>::impl
      { 
        // member fields
        const int n_dims;
        thrust_size_t n_part;
        detail::u01<real_t> rng;
        thrust_device::vector<real_t> 
          // particle attributes
          rd3, // dry radii cubed 
          xi, 
          x, 
          y, 
          z,
          // helper vectors
          u01 // uniform random numbers between 0 and 1
          ;

        // fills u01[0:n] with random numbers
        void urand(thrust_size_t n) { rng.generate_n(u01, n); }

        // compile-time min(1, n) 
        int m1(int n) { return n == 0 ? 1 : n; }

        // ctor 
        impl(real_t sd_conc_mean, int nx, int ny, int nz) : 
          n_dims(nx/m1(nx) + ny/m1(ny) + nz/m1(nz)), // 0, 1, 2 or 3
          n_part(sd_conc_mean * m1(nx) * m1(ny) * m1(nz)) // sd_conc_mean * nx * ny * nz (with ni=min(ni,1))
        {
          rd3.resize(n_part);
          u01.resize(n_part);
        }
      };

      // ctor
      template <typename real_t, int device>
      particles<real_t, device>::particles(
        real_t sd_conc_mean,
        int nx,
        int ny,
        int nz
      )
        : pimpl(new impl(sd_conc_mean, nx, ny, nz)) 
      { }
    };
  };
};
