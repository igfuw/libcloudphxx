/// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include <thrust/host_vector.h>
#include <thrust/iterator/constant_iterator.h>

#include <boost/array.hpp>
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
      bool init_called, should_now_run_async, selected_before_counting;

      // member fields
      opts_init_t<real_t> opts_init; // a copy
      const int n_dims;
      const thrust_size_t n_cell; 
      thrust_size_t n_part,            // total number of SDs
                    n_part_old,        // total number of SDs before source
                    n_part_to_init;    // number of SDs to be initialized by source
      detail::rng<real_t, device> rng;
      detail::config<real_t> config;

      // pointer to collision kernel
      kernel_base<real_t, n_t> *p_kernel;
 
      // containters for all kernel types
      thrust_device::vector<kernel_golovin<real_t, n_t> > k_golovin;
      thrust_device::vector<kernel_geometric<real_t, n_t> > k_geometric;
      thrust_device::vector<kernel_long<real_t, n_t> > k_long;
      thrust_device::vector<kernel_geometric_with_efficiencies<real_t, n_t> > k_geometric_with_efficiencies;
      thrust_device::vector<kernel_geometric_with_multiplier<real_t, n_t> > k_geometric_with_multiplier;
      thrust_device::vector<kernel_onishi<real_t, n_t> > k_onishi;

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
	x,   // x spatial coordinate (for 1D, 2D and 3D)
	y,   // y spatial coordinate (for 3D)
	z,   // z spatial coordinate (for 2D and 3D)
        sstp_tmp_rv, // either rv_old or advection-caused change in water vapour mixing ratio
        sstp_tmp_th, // ditto for theta_d
        sstp_tmp_rh; // ditto for rho

      // dry radii distribution characteristics
      real_t log_rd_min, // logarithm of the lower bound of the distr
             log_rd_max, // logarithm of the upper bound of the distr
             multiplier; // multiplier calculated for the above values

      // terminal velocity (per particle)
      thrust_device::vector<real_t> vt; 
      // sea level term velocity according to Beard 1977, compute once
      thrust_device::vector<real_t> vt_0; 

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

      // Eulerian-Lagrangian interface vars
      thrust_device::vector<real_t> 
        rhod,    // dry air density
        th,      // potential temperature (dry)
        rv,      // water vapour mixing ratio
        sstp_tmp_chem_0, // ditto for trace gases
        sstp_tmp_chem_1, // ditto for trace gases
        sstp_tmp_chem_2, // ditto for trace gases
        sstp_tmp_chem_3, // ditto for trace gases
        sstp_tmp_chem_4, // ditto for trace gases
        sstp_tmp_chem_5, // ditto for trace gases
        courant_x, 
        courant_y, 
        courant_z;

      std::map<enum chem_species_t, thrust_device::vector<real_t> > ambient_chem;

      // map of the accumulated volume/volume/mass of water/dry/chem that fell out of the domain
      std::map<enum output_t, real_t> output_puddle;
  
      thrust_device::vector<real_t> 
        T,  // temperature [K]
        p,  // pressure [Pa]
        RH, // relative humisity (p_v / p_vs)
        eta;// dynamic viscosity 

      // sorting needed only for diagnostics and coalescence
      bool sorted;

      // true if coalescence timestep has to be reduced, accesible from both device and host code
      bool *increase_sstp_coal;
      // is it a pure const_multi run, i.e. no sd_conc
      bool pure_const_multi;

      // timestep counter
      n_t stp_ctr;

      // maps linear Lagrangian component indices into Eulerian component linear indices
      // the map key is the address of the Thrust vector
      std::map<
        const thrust_device::vector<real_t>*, 
        thrust::host_vector<int> 
      > l2e; 

      // chem stuff
      // TODO: consider changing the unit to AMU or alike (very small numbers!)
      std::vector<typename thrust_device::vector<real_t>::iterator >
        chem_bgn, chem_end; // indexed with enum chem_species_t
      thrust_device::vector<real_t> chem_rhs, chem_ante_rhs, chem_post_rhs;
      /* TODO:
        On May 9, 2012, at 7:44 PM, Karsten Ahnert wrote:
        > ... unfortunately the Rosenbrock method cannot be used with any other state type than ublas.matrix.
        > ... I think, the best steppers for stiff systems and thrust are the
        > runge_kutta_fehlberg78 or the bulirsch_stoer with a very high order. But
        > should benchmark both steppers and choose the faster one.
      */
      boost::numeric::odeint::runge_kutta4<
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
        tmp_device_real_part1,  
        tmp_device_real_part2,  
        tmp_device_real_part3,
        tmp_device_real_part4,
        tmp_device_real_cell,
        tmp_device_real_cell1,
	&u01;  // uniform random numbers between 0 and 1 // TODO: use the tmp array as rand argument?
      thrust_device::vector<unsigned int>
        tmp_device_n_part,
        &un; // uniform natural random numbers between 0 and max value of unsigned int
      thrust_device::vector<thrust_size_t>
        tmp_device_size_cell,
        tmp_device_size_part;

      // to simplify foreach calls
      const thrust::counting_iterator<thrust_size_t> zero;

      // -- distributed memory stuff --
      // TODO: move to a separate struct?

      const int mpi_rank,
                mpi_size;

      // boundary type (shared mem/distmem)
      std::pair<detail::bcond_t, detail::bcond_t> bcond;

      // number of particles to be copied left/right in distmem setup
      unsigned int lft_count, rgt_count;

      // nx in devices to the left of this one
      unsigned int n_x_bfr,
                   n_x_tot; // total number of cells in x in all devices

      // number of cells in devices to the left of this one
      thrust_size_t n_cell_bfr;

      const int halo_x, // number of cells in the halo for courant_x before first "real" cell, halo only in x
                halo_y, // number of cells in the halo for courant_y before first "real" cell, halo only in x
                halo_z; // number of cells in the halo for courant_z before first "real" cell, halo only in x


      // x0 of the process to the right
      real_t rgt_x0;

      // x1 of the process to the left
      real_t lft_x1;

      // in/out buffers for SDs copied from other GPUs
      thrust_device::vector<n_t> in_n_bfr, out_n_bfr;
      // TODO: real buffers could be replaced with tmp_device_real_part1/2 if sstp_cond>1
      thrust_device::vector<real_t> in_real_bfr, out_real_bfr;

      // ids of sds to be copied with distmem
      thrust_device::vector<thrust_size_t> &lft_id, &rgt_id;

      // real_t vectors copied in distributed memory case
      std::vector<thrust_device::vector<real_t>*> distmem_real_vctrs;


      // methods

      // fills u01[0:n] with random numbers
      void rand_u01(thrust_size_t n) { rng.generate_n(u01, n); }

      // fills un[0:n] with random numbers
      void rand_un(thrust_size_t n) { rng.generate_n(un, n); }

      // max(1, n)
      int m1(int n) { return n == 0 ? 1 : n; }

      // ctor 
      impl(const opts_init_t<real_t> &_opts_init, const std::pair<detail::bcond_t, detail::bcond_t> &bcond, const int &mpi_rank, const int &mpi_size, const int &n_x_tot) : 
        init_called(false),
        should_now_run_async(false),
        selected_before_counting(false),
	opts_init(_opts_init),
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
        zero(0),
        n_part(0),
        sorted(false), 
        u01(tmp_device_real_part),
        n_user_params(opts_init.kernel_parameters.size()),
        un(tmp_device_n_part),
        rng(opts_init.rng_seed),
        stp_ctr(0),
	bcond(bcond),
        n_x_bfr(0),
        n_cell_bfr(0),
        mpi_rank(mpi_rank),
        mpi_size(mpi_size),
        lft_x1(-1),  // default to no
        rgt_x0(-1),  // MPI boudanry
        lft_id(i),   // note: reuses i vector
        rgt_id(tmp_device_size_part),
        n_x_tot(n_x_tot),
        halo_x( 
          n_dims == 1 ? 1:                 // 1D
            n_dims == 2 ? opts_init.nz:    // 2D
              opts_init.nz * opts_init.ny // 3D
        ),
        halo_y((opts_init.ny + 1) * opts_init.nz), // 3D
        halo_z( 
          n_dims == 2 ? opts_init.nz + 1:      // 2D
            (opts_init.nz + 1) * opts_init.ny),// 3D
        pure_const_multi (((opts_init.sd_conc) == 0) && (opts_init.sd_const_multi > 0 || opts_init.sd_const_multi_dry_sizes > 0)) // coal prob can be greater than one only in sd_conc simulations
      {

        // set 0 dev_count to mark that its not a multi_CUDA spawn
        // if its a spawn, multi_CUDA ctor will alter it
        opts_init.dev_count = 0; 

        // if using nvcc, put increase_sstp_coal flag in host memory, but with direct access from device code
#if defined(__NVCC__)
        cudaMallocHost(&increase_sstp_coal, sizeof(bool));
#else
        increase_sstp_coal = new bool();
#endif
        *increase_sstp_coal = false;

        // initialising host temporary arrays
        {
          thrust_size_t n_grid;
          switch (n_dims) // TODO: document that 3D is xyz, 2D is xz, 1D is x
          {
            case 3:
              n_grid = std::max(std::max(
                (opts_init.nx+2+1) * (opts_init.ny+0) * (opts_init.nz+0), 
                (opts_init.nx+2) * (opts_init.ny+1) * (opts_init.nz+0)),
                (opts_init.nx+2) * (opts_init.ny+0) * (opts_init.nz+1)
              );
              break;
            case 2:
              n_grid = std::max(
                (opts_init.nx+2+1) * (opts_init.nz+0), 
                (opts_init.nx+2) * (opts_init.nz+1)
              );
              break;
            case 1:
              n_grid = opts_init.nx+2+1;
              break;
            case 0:
              n_grid = 1;
              break;
            default: assert(false); 
          }
          if (n_dims != 0) assert(n_grid > n_cell);
	  tmp_host_real_grid.resize(n_grid);
        }

        typedef thrust_device::vector<real_t>* ptr_t;
        ptr_t arr[] = {&rd3, &rw2, &kpa, &vt, &x, &z};
        distmem_real_vctrs = std::vector<ptr_t>(arr, arr + sizeof(arr) / sizeof(ptr_t) );
        if(opts_init.ny > 0) distmem_real_vctrs.push_back(&y);
        if(opts_init.sstp_cond > 1 && opts_init.exact_sstp_cond)
        {
           distmem_real_vctrs.push_back(&sstp_tmp_rv);
           distmem_real_vctrs.push_back(&sstp_tmp_th);
           distmem_real_vctrs.push_back(&sstp_tmp_rh);
        }
      }

      void sanity_checks();
      void init_SD_with_distros();
      void init_SD_with_distros_sd_conc(const common::unary_function<real_t> &, const real_t &);
      void init_SD_with_distros_tail(const common::unary_function<real_t> &, const real_t);
      void init_SD_with_distros_const_multi(const common::unary_function<real_t> &);
      void init_SD_with_distros_finalize(const real_t &);
      void init_SD_with_sizes();
      void init_sanity_check(
        const arrinfo_t<real_t>, const arrinfo_t<real_t>, const arrinfo_t<real_t>,
        const arrinfo_t<real_t>, const arrinfo_t<real_t>, const arrinfo_t<real_t>,
        const std::map<enum chem_species_t, const arrinfo_t<real_t> >
      );

      void init_dry_sd_conc();
      void init_dry_const_multi(
        const common::unary_function<real_t> &n_of_lnrd
      );
      void init_dry_dry_sizes(real_t);

      void init_n_sd_conc(
        const common::unary_function<real_t> &n_of_lnrd
      );
      void init_n_const_multi(const thrust_size_t &);

      void dist_analysis_sd_conc(
        const common::unary_function<real_t> &n_of_lnrd,
        const n_t sd_conc,
        const real_t dt = 1.
      );
      void dist_analysis_const_multi(
        const common::unary_function<real_t> &n_of_lnrd 
      );
      void init_ijk();
      void init_xyz();
      void init_kappa(const real_t &);
      void init_count_num_sd_conc(const real_t & = 1);
      void init_count_num_const_multi(const common::unary_function<real_t> &);
      void init_count_num_const_multi(const common::unary_function<real_t> &, const thrust_size_t &);
      void init_count_num_dry_sizes(const real_t &);
      void init_count_num_hlpr(const real_t &, const thrust_size_t &);
      void init_e2l(const arrinfo_t<real_t> &, thrust_device::vector<real_t>*, const int = 0, const int = 0, const int = 0, const long int = 0);
      void init_wet();
      void init_sync();
      void init_grid();
      void init_hskpng_ncell();
      void init_hskpng_npart();
      void init_chem();
      void init_chem_aq();
      void init_sstp();
      void init_sstp_chem();
      void init_kernel();
      void init_vterm();

      void fill_outbuf();
      void mpi_exchange();

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
      void hskpng_resize_npart();

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
        const real_t power,
        const bool specific = true
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
      template<class adve_t>
      void adve_calc(bool, thrust_size_t = 0);
      void sedi();

      void cond_dm3_helper();
      void cond(const real_t &dt, const real_t &RH_max);
      void cond_sstp(const real_t &dt, const real_t &RH_max);
      void update_th_rv(thrust_device::vector<real_t> &);
      void update_state(thrust_device::vector<real_t> &, thrust_device::vector<real_t> &);
      void update_pstate(thrust_device::vector<real_t> &, thrust_device::vector<real_t> &);

      void coal(const real_t &dt);

      void chem_vol_ante();
      void chem_flag_ante();
      void chem_henry(const real_t &dt);
      void chem_dissoc();
      void chem_react(const real_t &dt);
      void chem_cleanup();
 
      thrust_size_t rcyc();
      void bcnd();

      void src(const real_t &dt);

      void sstp_step(const int &step, const bool &var_rho);
      void sstp_step_exact(const int &step, const bool &var_rho);
      void sstp_save();
      void sstp_step_chem(const int &step, const bool &var_rho);
      void sstp_save_chem();

      void post_copy(const opts_t<real_t>&);

      // distmem stuff
      void xchng_domains();
      bool distmem_mpi();
      bool distmem_cuda();
      bool distmem();
      void pack_n_lft();
      void pack_n_rgt();
      void pack_real_lft();
      void pack_real_rgt();
      void unpack_n(const int &);
      void unpack_real(const int &);
      void flag_lft();
      void flag_rgt();
      void bcnd_remote_lft(const real_t &, const real_t &);
      void bcnd_remote_rgt(const real_t &, const real_t &);
    };
  };
};
