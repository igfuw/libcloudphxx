// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    // pimpl stuff
    template <typename real_t, int device>
    struct particles<real_t, device>::impl
    { 
      typedef unsigned long n_t;

      // member fields
      const int n_dims;
      const thrust_size_t n_part; 
      const real_t sd_conc_mean; //, dx, dy, dz;
      detail::u01<real_t, device> rng;
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
      thrust_device::vector<n_t>
	n;  // multiplicity

      // fills u01[0:n] with random numbers
      void urand(thrust_size_t n) { rng.generate_n(u01, n); }

      // compile-time min(1, n) 
      int m1(int n) { return n == 0 ? 1 : n; }

      // ctor 
      impl(
	real_t sd_conc_mean, 
	int nx, int ny, int nz
      ) : 
	sd_conc_mean(sd_conc_mean),
	n_dims(nx/m1(nx) + ny/m1(ny) + nz/m1(nz)), // 0, 1, 2 or 3
	n_part(sd_conc_mean * m1(nx) * m1(ny) * m1(nz)) // sd_conc_mean * nx * ny * nz (with ni=min(ni,1))
      {
std::cerr << "device = " << device << std::endl;
#define str_value(arg) #arg
#define name(arg) str_value(arg)
#define thrust_device_name name(thrust_device)
std::cerr << "thrust_device = " << std::string(thrust_device_name) << std::endl;
	u01.resize(n_part);
	rd3.resize(n_part);
	n.resize(n_part);
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
