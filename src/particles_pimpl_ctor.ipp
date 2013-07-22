// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include <thrust/host_vector.h>

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
      const opts_t<real_t> opts; // nx, ny, nz, dx, dy, dz, ...;
      const int n_dims;
      const int n_cell; 
      const thrust_size_t n_part; 
      detail::u01<real_t, device> rng;

      // particle attributes
      thrust_device::vector<real_t> 
	rd3, // dry radii cubed 
	xi,  // wet radius proxy variable (e.g. xi = ln(rw))
        kpa, // kappa
	x,   // x spatial coordinate (for 2D and 3D)
	y,   // y spatial coordinate (for 3D)
	z;   // z spatial coordinate (for 1D, 2D and 3D)
      thrust_device::vector<n_t>
	n;  // multiplicity

      // helper vectors
      thrust_device::vector<real_t>
	u01; // uniform random numbers between 0 and 1

      // Eulerian-Lagrangian interface vers
      thrust_device::vector<int> 
        i, j, k, ijk; // Eulerian grid cell indices (always zero for 0D)
      thrust_device::vector<real_t> 
        rhod,    // dry air density
        rhod_th, // energy volume density divided by c_p
        rhod_rv, // water vapour density (=rhod * r_v)
        T, // temperature [K]
        p, // pressure [Pa]
        r; // water vapour mixing ratio [kg/kg]

      thrust::host_vector<thrust_size_t> 
        l2e; // maps linear Lagrangian component indices into Eulerian component linear indices
      thrust::host_vector<real_t>
        tmp_host_real_cell;

      // to simplify foreach calls
      const thrust::counting_iterator<thrust_size_t> zero;

      // fills u01[0:n] with random numbers
      void urand(thrust_size_t n) { rng.generate_n(u01, n); }

      // compile-time min(1, n) 
      int m1(int n) { return n == 0 ? 1 : n; }

      // ctor 
      impl(const opts_t<real_t> opts) : 
	opts(opts),
	n_dims(opts.nx/m1(opts.nx) + opts.ny/m1(opts.ny) + opts.nz/m1(opts.nz)), // 0, 1, 2 or 3
        n_cell(m1(opts.nx) * m1(opts.ny) * m1(opts.nz)),
	n_part(opts.sd_conc_mean * n_cell), // TODO: what if multiple spectra/kappas
        zero(0)
      {
	u01.resize(n_part);
        i.resize(n_part); // 
        j.resize(n_part); //  > TODO: are they needed at all?
        k.resize(n_part); // 
        ijk.resize(n_part);
        tmp_host_real_cell.resize(n_cell);
      }

      // methods
      void sanity_checks();

      void init_dry(const unary_function<real_t> *n_of_lnrd);
      void init_xyz();
      void init_e2l(const ptrdiff_t *);
      void init_Tpr();
      void init_wet();

      void hskpng_ijk();
      void hskpng_Tpr();
      void sync(
        real_t *, // from // TODO: const
        thrust_device::vector<real_t> & // to
      );
      void sync(
        const thrust_device::vector<real_t> &, // from
        real_t * // to
      );
    };

    // ctor
    template <typename real_t, int device>
    particles<real_t, device>::particles(const opts_t<real_t> opts) :
      pimpl(new impl(opts))
    {
      pimpl->sanity_checks();
      assert(opts.dx != 0);
      assert(opts.dy != 0);
      assert(opts.dz != 0);
    }
  };
};
