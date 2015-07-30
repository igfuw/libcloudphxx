// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include <thrust/host_vector.h>
#include <thrust/iterator/constant_iterator.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>

#include <map>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      enum { invalid = -1 };
    };  

    // pimpl stuff 
    template <typename real_t, backend_t device>
    struct particles_t<real_t, device>::impl
    { 
      // CUDA does not support max(unsigned long, unsigned long) -> using unsigned long long
      typedef unsigned long long n_t; // thrust_size_t?
 
      // order of operation flags
      bool should_now_run_async, selected_before_counting;

      // member fields
      const opts_init_t<real_t> opts_init; // a copy
      const int n_dims;
      const int n_cell; 
      thrust_size_t n_part; 
      detail::rng<real_t, device> rng;

      // pointer to collision kernel
      kernel_base<real_t, n_t> *p_kernel;
 
      //containters for all kernel types
      thrust_device::vector<kernel_golovin<real_t, n_t> > k_golovin;
      thrust_device::vector<kernel_geometric<real_t, n_t> > k_geometric;
      thrust_device::vector<kernel_long<real_t, n_t> > k_long;
      thrust_device::vector<kernel_geometric_with_efficiencies<real_t, n_t> > k_geometric_with_efficiencies;
      thrust_device::vector<kernel_geometric_with_multiplier<real_t, n_t> > k_geometric_with_multiplier;

      // device container for kernel parameters, could come from opts_init or a file depending on the kernel
      thrust_device::vector<real_t> kernel_parameters;

      //number of kernel parameters defined by user in opts_init
      const n_t n_user_params;

      // particle attributes
      thrust_device::vector<n_t>
	n;   // multiplicity
      thrust_device::vector<real_t> 
	rd3, // dry radii cubed 
	rw2, // wet radius square
        kpa, // kappa
	x,   // x spatial coordinate (for 2D and 3D)
	y,   // y spatial coordinate (for 3D)
	z;   // z spatial coordinate (for 1D, 2D and 3D)

      // terminal velocity (per particle)
      thrust_device::vector<real_t> vt; 

      // grid-cell volumes (per grid cell)
      thrust_device::vector<real_t> dv;

      // housekeeping data (per particle)
      thrust_device::vector<thrust_size_t> 
        i, j, k, ijk, // Eulerian grid cell indices (always zero for 0D)
        sorted_id, sorted_ijk;

      // Arakawa-C grid helper vars
      thrust_device::vector<thrust_size_t> 
        lft, rgt, abv, blw, fre, hnd; // TODO: could be reused after advection!

      // moment-counting stuff
      thrust_device::vector<thrust_size_t> 
        count_ijk; // key-value pair for sorting particles by cell index
      thrust_device::vector<n_t>
        count_num; // number of particles in a given grid cell
      thrust_device::vector<real_t> 
        count_mom; // statistical moment // TODO (perhaps tmp_device_real_cell could be referenced?)
      thrust_size_t count_n;

      // Eulerian-Lagrangian interface vers
      thrust_device::vector<real_t> 
        rhod,    // dry air density
        th,      // potential temperature (dry)
        rv,      // water vapour mixing ratio
        sstp_tmp_rv, // either rv_old or advection-caused change in water vapour mixing ratio
        sstp_tmp_th, // ditto for theta_d
        sstp_tmp_rh, // ditto for rho
        courant_x, 
        courant_y, 
        courant_z;
  
      thrust_device::vector<real_t> 
        T,  // temperature [K]
        p,  // pressure [Pa]
        RH, // relative humisity (p_v / p_vs)
        eta;// dynamic viscosity 

      // sorting needed only for diagnostics and coalescence
      bool sorted;

      // maps linear Lagrangian component indices into Eulerian component linear indices
      // the map key is the address of the Thrust vector
      std::map<
        const thrust_device::vector<real_t>*, 
        thrust::host_vector<thrust_size_t> 
      > l2e; 

      // chem stuff
      // TODO: consider changing the unit to AMU or alike (very small numbers!)
      std::vector<typename thrust_device::vector<real_t>::iterator >
        chem_bgn, chem_end; // indexed with enum chem_species_t
      thrust_device::vector<real_t> chem_noneq, chem_equil;
      /* TODO:
        On May 9, 2012, at 7:44 PM, Karsten Ahnert wrote:
        > ... unfortunately the Rosenbrock method cannot be used with any other state type than ublas.matrix.
        > ... I think, the best steppers for stiff systems and thrust are the
        > runge_kutta_fehlberg78 or the bulirsch_stoer with a very high order. But
        > should benchmark both steppers and choose the faster one.
      */
      boost::numeric::odeint::euler<
        thrust_device::vector<real_t>, // state_type
        real_t,                        // value_type
        thrust_device::vector<real_t>, // deriv_type
        real_t,                        // time_type
        boost::numeric::odeint::thrust_algebra,
        boost::numeric::odeint::thrust_operations,
        boost::numeric::odeint::never_resizer
      > chem_stepper;

      // temporary data
      thrust::host_vector<real_t>
        tmp_host_real_grid,
        tmp_host_real_cell;
      thrust::host_vector<thrust_size_t>
        tmp_host_size_cell;
      thrust_device::vector<real_t>
        tmp_device_real_part,
        tmp_device_real_cell,
	&u01;  // uniform random numbers between 0 and 1 // TODO: use the tmp array as rand argument?
      thrust_device::vector<unsigned int>
        tmp_device_n_part,
        &un; // uniform natural random numbers between 0 and max value of unsigned int
      thrust_device::vector<thrust_size_t>
        tmp_device_size_cell;

      // to simplify foreach calls
      const thrust::counting_iterator<thrust_size_t> zero;

      // fills u01[0:n] with random numbers
      void rand_u01(thrust_size_t n) { rng.generate_n(u01, n); }

      // fills un[0:n] with random numbers
      void rand_un(thrust_size_t n) { rng.generate_n(un, n); }

      // compile-time min(1, n) 
      int m1(int n) { return n == 0 ? 1 : n; }

      // ctor 
      impl(const opts_init_t<real_t> &opts_init) : 
        should_now_run_async(false),
        selected_before_counting(false),
	opts_init(opts_init),
	n_dims( // 0, 1, 2 or 3
          opts_init.nx/m1(opts_init.nx) + 
          opts_init.ny/m1(opts_init.ny) + 
          opts_init.nz/m1(opts_init.nz)
        ), 
        n_cell(
          m1(opts_init.nx) * 
          m1(opts_init.ny) *
          m1(opts_init.nz)
        ),
	n_part( // TODO: what if multiple spectra/kappas
          opts_init.sd_conc_mean * 
	  ((opts_init.x1 - opts_init.x0) / opts_init.dx) *
	  ((opts_init.y1 - opts_init.y0) / opts_init.dy) *
	  ((opts_init.z1 - opts_init.z0) / opts_init.dz)
        ),
        zero(0), 
        sorted(false), 
        u01(tmp_device_real_part),
        n_user_params(opts_init.kernel_parameters.size()),
        un(tmp_device_n_part),
        rng(opts_init.rng_seed)
      {
        // sanity checks
        if (n_dims > 0)
        {
	  if (!(opts_init.x0 >= 0 && opts_init.x0 < m1(opts_init.nx) * opts_init.dx))
            throw std::runtime_error("!(x0 >= 0 & x0 < min(1,nx)*dz)"); 
	  if (!(opts_init.y0 >= 0 && opts_init.y0 < m1(opts_init.ny) * opts_init.dy))
            throw std::runtime_error("!(y0 >= 0 & y0 < min(1,ny)*dy)"); 
	  if (!(opts_init.z0 >= 0 && opts_init.z0 < m1(opts_init.nz) * opts_init.dz))
            throw std::runtime_error("!(z0 >= 0 & z0 < min(1,nz)*dz)"); 
	  if (!(opts_init.x1 > opts_init.x0 && opts_init.x1 <= m1(opts_init.nx) * opts_init.dx))
            throw std::runtime_error("!(x1 > x0 & x1 <= min(1,nx)*dx)");
	  if (!(opts_init.y1 > opts_init.y0 && opts_init.y1 <= m1(opts_init.ny) * opts_init.dy))
            throw std::runtime_error("!(y1 > y0 & y1 <= min(1,ny)*dy)");
	  if (!(opts_init.z1 > opts_init.z0 && opts_init.z1 <= m1(opts_init.nz) * opts_init.dz))
            throw std::runtime_error("!(z1 > z0 & z1 <= min(1,nz)*dz)");
        }

        if (opts_init.dt == 0) throw std::runtime_error("please specify opts_init.dt");
        if (opts_init.sd_conc_mean == 0) throw std::runtime_error("please specify opts_init.sd_conc");
        if (opts_init.coal_switch)
        {
          if(opts_init.terminal_velocity == vt_t::undefined) throw std::runtime_error("please specify opts_init.terminal_velocity or turn off opts_init.coal_switch");
          if(opts_init.kernel == kernel_t::undefined) throw std::runtime_error("please specify opts_init.kernel");
        }
        if (opts_init.sedi_switch)
          if(opts_init.terminal_velocity == vt_t::undefined) throw std::runtime_error("please specify opts_init.terminal_velocity or turn off opts_init.sedi_switch");

        // note: there could be less tmp data spaces if _cell vectors
        //       would point to _part vector data... but using.end() would not possible
        // initialising device temporary arrays
	tmp_device_real_part.resize(n_part);
        tmp_device_real_cell.resize(n_cell);
        tmp_device_size_cell.resize(n_cell);
	tmp_device_n_part.resize(n_part);

        // initialising host temporary arrays
        {
          int n_grid;
          switch (n_dims) // TODO: document that 3D is xyz, 2D is xz, 1D is z
          {
            case 3:
              n_grid = std::max(std::max(
                (opts_init.nx+1) * (opts_init.ny+0) * (opts_init.nz+0), 
                (opts_init.nx+0) * (opts_init.ny+1) * (opts_init.nz+0)),
                (opts_init.nx+0) * (opts_init.ny+0) * (opts_init.nz+1)
              );
              break;
            case 2:
              n_grid = std::max(
                (opts_init.nx+1) * (opts_init.nz+0), 
                (opts_init.nx+0) * (opts_init.nz+1)
              );
              break;
            case 0:
              n_grid = 1;
              break;
            default: assert(false); // TODO: 1D case
          }
          if (n_dims != 0) assert(n_grid > n_cell);
	  tmp_host_real_grid.resize(n_grid);
        }
        tmp_host_size_cell.resize(n_cell);
        tmp_host_real_cell.resize(n_cell);
      }

      // methods
      void sanity_checks();

      void init_dry(
        const real_t kappa, // TODO: map
        const common::unary_function<real_t> *n_of_lnrd
      );
      void init_xyz();
      void init_e2l(const arrinfo_t<real_t> &, thrust_device::vector<real_t>*, const int = 0, const int = 0, const int = 0);
      void init_wet();
      void init_sync();
      void init_grid();
      void init_hskpng();
      void init_chem();
      void init_sstp();
      void init_kernel();

      void fill_outbuf();

           // rename hskpng_ -> step_?
      void hskpng_sort_helper(bool);
      void hskpng_sort();
      void hskpng_shuffle_and_sort();
      void hskpng_count();
      void hskpng_ijk();
      void hskpng_Tpr();

      void hskpng_vterm_all();
      void hskpng_vterm_invalid();
      void hskpng_remove_n0();

      void moms_all();
   
      void moms_cmp(
        const typename thrust_device::vector<real_t>::iterator &vec1_bgn,
        const typename thrust_device::vector<real_t>::iterator &vec2_bgn
      );
      void moms_ge0(
        const typename thrust_device::vector<real_t>::iterator &vec_bgn
      );
      void moms_rng(
        const real_t &min, const real_t &max, 
        const typename thrust_device::vector<real_t>::iterator &vec_bgn
      ); 
      void moms_calc(
	const typename thrust_device::vector<real_t>::iterator &vec_bgn,
        const real_t power
      );

      void mass_dens_estim(
	const typename thrust_device::vector<real_t>::iterator &vec_bgn,
        const real_t, const real_t, const real_t
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
      void sedi();

      void cond_dm3_helper();
      void cond(const real_t &dt, const real_t &RH_max);

      void coal(const real_t &dt);

      void chem(const real_t &dt, const std::vector<real_t> &chem_gas, 
                const bool &chem_dsl, const bool &chem_dsc, const bool &chem_rct);
      thrust_size_t rcyc();
      real_t bcnd(); // returns accumulated rainfall

      void sstp_step(const int &step, const bool &var_rho);
      void sstp_save();
    };

    // ctor
    template <typename real_t, backend_t device>
    particles_t<real_t, device>::particles_t(const opts_init_t<real_t> &opts_init) :
      pimpl(new impl(opts_init))
    {
      this->opts_init = &pimpl->opts_init;
      pimpl->sanity_checks();
    }

    // outbuf
    template <typename real_t, backend_t device>
    real_t *particles_t<real_t, device>::outbuf() 
    {
      pimpl->fill_outbuf();
      return &(*(pimpl->tmp_host_real_cell.begin()));
    }
  };
};
