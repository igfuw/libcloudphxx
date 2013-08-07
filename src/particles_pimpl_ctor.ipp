// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include <thrust/host_vector.h>
#include <map>

namespace libcloudphxx
{
  namespace lgrngn
  {
 
    // pimpl stuff 
    template <typename real_t, int device>
    struct particles<real_t, device>::impl
    { 
      typedef unsigned long n_t; // thrust_size_t?

      // member fields
      const opts_t<real_t> opts; // a copy
      const int n_dims;
      const int n_cell; 
      const thrust_size_t n_part; 
      detail::u01<real_t, device> rng;

      // particle attributes
      thrust_device::vector<real_t> 
	rd3, // dry radii cubed 
	rw2, // wet radius square
        kpa, // kappa
	x,   // x spatial coordinate (for 2D and 3D)
	y,   // y spatial coordinate (for 3D)
	z;   // z spatial coordinate (for 1D, 2D and 3D)
      thrust_device::vector<n_t>
	n;  // multiplicity

      //
      thrust_device::vector<thrust_size_t> 
        i, j, k, ijk, // Eulerian grid cell indices (always zero for 0D)
        sorted_id, sorted_ijk;

      // 2D Arakawa-C grid helper vars
      thrust_device::vector<thrust_size_t> 
        lft, rgt, abv, blw;

      //
      thrust_device::vector<thrust_size_t> 
        count_ijk; // key-value pair for sorting particles by cell index
      thrust_device::vector<n_t>
        count_num; // number of particles in a given grid cell
      thrust_device::vector<real_t> 
        count_mom; // statistical moment
      thrust_size_t count_n;

      // Eulerian-Lagrangian interface vers
      thrust_device::vector<real_t> 
        rhod,    // dry air density
        rhod_th, // energy volume density divided by c_p
        rhod_rv, // water vapour density (=rhod * r_v)
        courant_x, 
        courant_y, 
        courant_z;
  
      thrust_device::vector<real_t> 
        T, // temperature [K]
        p, // pressure [Pa]
        r; // water vapour mixing ratio [kg/kg]

      // sorting needed only for diagnostics and coalescence
      bool sorted;


      // maps linear Lagrangian component indices into Eulerian component linear indices
      // the map key is the address of the Thrust vector
      std::map<
        const thrust_device::vector<real_t>*, 
        thrust::host_vector<thrust_size_t> 
      > l2e; 

      // temporary data
      thrust::host_vector<real_t>
        tmp_host_real_grid,
        &tmp_host_real_cell;
      thrust::host_vector<thrust_size_t>
        tmp_host_size_cell;

      // helper vectors
      thrust_device::vector<real_t>
        tmp_device_real_part,
	&u01; // uniform random numbers between 0 and 1 // TODO: use the tmp array as rand argument?

      // to simplify foreach calls
      const thrust::counting_iterator<thrust_size_t> zero;

      // fills u01[0:n] with random numbers
      void rand_u01(thrust_size_t n) { rng.generate_n(u01, n); }

      // compile-time min(1, n) 
      int m1(int n) { return n == 0 ? 1 : n; }

      // ctor 
      impl(const opts_t<real_t> &opts) : 
	opts(opts),
	n_dims(opts.nx/m1(opts.nx) + opts.ny/m1(opts.ny) + opts.nz/m1(opts.nz)), // 0, 1, 2 or 3
        n_cell(m1(opts.nx) * m1(opts.ny) * m1(opts.nz)),
	n_part(opts.sd_conc_mean * n_cell), // TODO: what if multiple spectra/kappas
        zero(0), // TODO: is it used anywhere?
        sorted(false), 
        tmp_host_real_cell(tmp_host_real_grid),
        u01(tmp_device_real_part)
      {
        // initialising device temporary arrays
	tmp_device_real_part.resize(n_part);

        // initialising host temporary arrays
        {
          int n_grid;
          switch (n_dims)
          {
            case 0: 
              break;
            case 2:
              n_grid = std::max((opts.nx+1)*opts.nz, opts.nx*(opts.nz+1));
              break;
            default: assert(false);
          }
          assert(n_grid > n_cell);
	  tmp_host_real_grid.resize(n_grid);
        }
        tmp_host_size_cell.resize(n_cell);
      }

      // methods
      void sanity_checks();

      void init_dry(const common::unary_function<real_t> *n_of_lnrd);
      void init_xyz();
      void init_e2l(const arrinfo_t<real_t> &, thrust_device::vector<real_t>*, const int = 0, const int = 0, const int = 0);
      void init_wet();
      void init_sync();
      void init_grid();
      void init_hskpng();

      void fill_outbuf();

           // rename hskpng_ -> step_?
      void hskpng_shuffle();
      void hskpng_sort();
      void hskpng_count();
      void hskpng_ijk();
      void hskpng_Tpr();

      void moms_rng(
        const real_t &min, const real_t &max, 
        const thrust_device::vector<real_t> &radii
      ); 
      void moms_calc(
	const thrust_device::vector<real_t> &radii,
	const real_t power
      );

      void sync(
        const arrinfo_t<real_t> &, // from 
        thrust_device::vector<real_t> & // to
      );
      void sync(
        const thrust_device::vector<real_t> &, // from
        arrinfo_t<real_t> &// to
      );

      void adve();
    };

    // ctor
    template <typename real_t, int device>
    particles<real_t, device>::particles(const opts_t<real_t> &opts) :
      pimpl(new impl(opts))
    {
      pimpl->sanity_checks();
    }

    // outbuf
    template <typename real_t, int device>
    real_t *particles<real_t, device>::outbuf() 
    {
      pimpl->fill_outbuf();
      return &(*(pimpl->tmp_host_real_cell.begin()));
    }
  };
};
