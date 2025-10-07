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
//#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
//#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include <map>
#include <set>


namespace libcloudphxx
{
  namespace lgrngn
  {
    // pimpl stuff 
    template <typename real_t, backend_t device>
    struct particles_t<real_t, device>::impl
    { 
      // CUDA does not support max(unsigned long, unsigned long) -> using unsigned long long
      typedef unsigned long long n_t; // thrust_size_t?
 
      // order of operation flags
      bool init_called, should_now_run_async, selected_before_counting, should_now_run_cond;

      // did density vary in this step
      bool var_rho;

      // member fields
      opts_init_t<real_t> opts_init; // a copy
      const int n_dims;
      const thrust_size_t n_cell; 
      thrust_size_t n_part,            // total number of SDs
                    n_part_old,        // total number of SDs before source
                    n_part_to_init;    // number of SDs to be initialized by source
      detail::rng<real_t, device> rng;
      detail::config<real_t> config;
      as_t adve_scheme;         // actual advection scheme used, might be different from opts_init.adve_scheme if courant>halo

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
        up,  // turbulent perturbation of velocity
        vp,  // turbulent perturbation of velocity
        wp,  // turbulent perturbation of velocity
        ssp, // turbulent perturbation of supersaturation
        dot_ssp, // time derivative of the turbulent perturbation of supersaturation
        sstp_tmp_rv, // either rv_old or advection-caused change in water vapour mixing ratio
        sstp_tmp_th, // ditto for theta
        sstp_tmp_rh, // ditto for rho
        sstp_tmp_p, // ditto for pressure
        incloud_time; // time this SD has been within a cloud

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
      std::map<enum common::output_t, real_t> output_puddle;
  
      thrust_device::vector<real_t> 
        T,  // temperature [K]
        p,  // pressure [Pa]
        RH, // relative humisity 
        eta,// dynamic viscosity 
        diss_rate; // turbulent kinetic energy dissipation rate

      thrust_device::vector<real_t> w_LS; // large-scale subsidence velocity profile
      thrust_device::vector<real_t> SGS_mix_len; // SGS mixing length profile
      thrust_device::vector<real_t> aerosol_conc_factor; // profile of aerosol concentration factor

      // time steps to be used, considering that opts.dt can override opts_init.dt
      real_t dt;
      int sstp_cond, sstp_coal, sstp_chem;

      // sorting needed only for diagnostics and coalescence
      bool sorted;

      // true if coalescence timestep has to be reduced, accesible from both device and host code
      bool *increase_sstp_coal;
      // is it a pure const_multi run, i.e. no sd_conc
      bool pure_const_multi;

      // is it allowed to do substepping, if not, some memory can be saved
      bool allow_sstp_cond,
           allow_sstp_chem;

      // timestep counter
      n_t src_stp_ctr, rlx_stp_ctr;

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
      tmp_vector_pool<thrust::host_vector<real_t>> 
        tmp_host_real_grid,
        tmp_host_real_cell;
      tmp_vector_pool<thrust::host_vector<thrust_size_t>> 
        tmp_host_size_cell;
      tmp_vector_pool<thrust_device::vector<real_t>>       
        tmp_device_real_part,
        // tmp_device_real_part1,  
        // tmp_device_real_part2,  
        // tmp_device_real_part3,
        // tmp_device_real_part4,
        // tmp_device_real_part5,
        tmp_device_real_cell;
        // tmp_device_real_cell1,
        // tmp_device_real_cell2,
      tmp_vector_pool<thrust_device::vector<unsigned int>>
        tmp_device_n_part;
      tmp_vector_pool<thrust_device::vector<thrust_size_t>>
        tmp_device_size_cell,
        tmp_device_size_part;

      // guards for temp vectors that are used in multiple functions and need to stay unchanged inbetween
      std::uniqe_ptr<
        tmp_vector_pool<thrust_device::vector<real_t>>::guard
      > n_filtered_gp,
        V_gp;


      // to simplify foreach calls
      const thrust::counting_iterator<thrust_size_t> zero;

      // -- distributed memory stuff --
      // TODO: move to a separate struct?

      const int mpi_rank,
                mpi_size;

      // boundary type in x direction (shared mem/distmem/open/periodic)
      std::pair<detail::bcond_t, detail::bcond_t> bcond;

      // number of particles to be copied left/right in distmem setup
      unsigned int lft_count, rgt_count;

      // nx in devices to the left of this one
      unsigned int n_x_bfr,
                   n_x_tot; // total number of cells in x in all devices of this process

      // number of cells in devices to the left of this one
      thrust_size_t n_cell_bfr;

      const int halo_size, // NOTE:  halo_size = 0 means that both x courant numbers in edge cells are known, what is equivalent to halo = 1 in libmpdata++
                           // NOTE2: halo means that some values of the Eulerian courant array are pointed to by more than one e2l, what could lead to race conditions if we wanted to sync out courants
                halo_x, // number of cells in the halo for courant_x before first "real" cell, halo only in x
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

      // --- containters with vector pointers to help resize and copy vectors ---

      // vectors copied between distributed memories (MPI, multi_CUDA), these are SD attributes
      std::set<std::pair<thrust_device::vector<real_t>*, real_t>>         distmem_real_vctrs; // pair of vector and its initial value
      std::set<thrust_device::vector<n_t>*>                               distmem_n_vctrs;
//      std::set<thrust_device::vector<thrust_size_t>*>  distmem_size_vctrs; // no size vectors copied?
//
      // vetors that are not in distmem_real_vctrs that need to be resized when the number of SDs changes, these are helper variables
//      std::set<thrust_device::vector<real_t>*>         resize_real_vctrs;
//      std::set<thrust_device::vector<n_t>*>            resize_n_vctrs;
      std::set<thrust_device::vector<thrust_size_t>*>  resize_size_vctrs;


      // --- methods ---

      // fills u01 with n random real numbers uniformly distributed in range [0,1)
      void rand_u01(thrust_device::vector<real_t> &u01, thrust_size_t n) { rng.generate_n(u01, n); }

      // fills un with n random integers uniformly distributed on the whole integer range
      void rand_un(thrust_device::vector<unsigned int> &un, thrust_size_t n) { rng.generate_n(un, n); }

      // max(1, n)
      int m1(int n) { return n == 0 ? 1 : n; }

      // ctor 
      impl(const opts_init_t<real_t> &_opts_init, const std::pair<detail::bcond_t, detail::bcond_t> &bcond, const int &mpi_rank, const int &mpi_size, const int &n_x_tot) : 
        init_called(false),
        should_now_run_async(false),
        selected_before_counting(false),
        should_now_run_cond(false),
        var_rho(false),
        opts_init(_opts_init),
        n_dims( // 0, 1, 2 or 3
          _opts_init.nx/m1(_opts_init.nx) + 
          _opts_init.ny/m1(_opts_init.ny) + 
          _opts_init.nz/m1(_opts_init.nz)
        ), 
        n_cell(
          m1(_opts_init.nx) * 
          m1(_opts_init.ny) *
          m1(_opts_init.nz)
        ),
        zero(0),
        n_part(0),
        sorted(false), 
        n_user_params(_opts_init.kernel_parameters.size()),
        rng(_opts_init.rng_seed),
        src_stp_ctr(0),
        rlx_stp_ctr(0),
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
        halo_size(_opts_init.adve_scheme == as_t::pred_corr ? 2 : 0), 
        halo_x( 
          n_dims == 1 ? halo_size:                                      // 1D
          n_dims == 2 ? halo_size * _opts_init.nz:                       // 2D
                        halo_size * _opts_init.nz * _opts_init.ny         // 3D
        ),
        halo_y(         halo_size * (_opts_init.ny + 1) * _opts_init.nz), // 3D
        halo_z( 
          n_dims == 2 ? halo_size * (_opts_init.nz + 1):                 // 2D
                        halo_size * (_opts_init.nz + 1) * _opts_init.ny   // 3D
        ),
        w_LS(_opts_init.w_LS),
        SGS_mix_len(_opts_init.SGS_mix_len),
        aerosol_conc_factor(_opts_init.aerosol_conc_factor),
        adve_scheme(_opts_init.adve_scheme),
        allow_sstp_cond(_opts_init.sstp_cond > 1 || _opts_init.variable_dt_switch),
        allow_sstp_chem(_opts_init.sstp_chem > 1 || _opts_init.variable_dt_switch),
        pure_const_multi (((_opts_init.sd_conc) == 0) && (_opts_init.sd_const_multi > 0 || _opts_init.dry_sizes.size() > 0)), // coal prob can be greater than one only in sd_conc simulations
        //tmp_device_real_part(6),
        tmp_device_real_cell(3) // 3 temporary vectors of this type; NOTE: default is 1
      {

        // set 0 dev_count to mark that its not a multi_CUDA spawn
        // if its a spawn, multi_CUDA ctor will alter it
        opts_init.dev_count = 0; 

        // if using nvcc, put increase_sstp_coal flag in host memory, but with direct access from device code
#if defined(__NVCC__)
        gpuErrchk(cudaMallocHost(&increase_sstp_coal, sizeof(bool)));
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
                (opts_init.nx+2*halo_size+1) * (opts_init.ny+0) * (opts_init.nz+0), 
                (opts_init.nx+2*halo_size) * (opts_init.ny+1) * (opts_init.nz+0)),
                (opts_init.nx+2*halo_size) * (opts_init.ny+0) * (opts_init.nz+1)
              );
              break;
            case 2:
              n_grid = std::max(
                (opts_init.nx+2*halo_size+1) * (opts_init.nz+0), 
                (opts_init.nx+2*halo_size) * (opts_init.nz+1)
              );
              break;
            case 1:
              n_grid = opts_init.nx+2*halo_size+1;
              break;
            case 0:
              n_grid = 1;
              break;
            default: assert(false); 
          }
          if (n_dims != 0) assert(n_grid > n_cell);
          tmp_host_real_grid.resize(n_grid);
        }

        // initializing distmem_real_vctrs - list of real_t vectors with properties of SDs that have to be copied/removed/recycled when a SD is copied/removed/recycled
        // NOTE: this does not include chemical stuff due to the way chem vctrs are organized! multi_CUDA / MPI does not work with chemistry as of now
        distmem_real_vctrs.insert({&rd3, detail::no_initial_value});
        distmem_real_vctrs.insert({&rw2, detail::no_initial_value});
        distmem_real_vctrs.insert({&kpa, detail::no_initial_value});

        distmem_real_vctrs.insert({&vt,  detail::invalid});

        if (opts_init.nx != 0)  distmem_real_vctrs.insert({&x, detail::no_initial_value});
        if (opts_init.ny != 0)  distmem_real_vctrs.insert({&y, detail::no_initial_value});
        if (opts_init.nz != 0)  distmem_real_vctrs.insert({&z, detail::no_initial_value});

        if(allow_sstp_cond && opts_init.exact_sstp_cond)
        {
           distmem_real_vctrs.insert({&sstp_tmp_rv, detail::no_initial_value});
           distmem_real_vctrs.insert({&sstp_tmp_th, detail::no_initial_value});
           distmem_real_vctrs.insert({&sstp_tmp_rh, detail::no_initial_value});
           // sstp_tmp_p needs to be added if a constant pressure profile is used, but this is only known after init - see particles_init
        }

        if(opts_init.turb_adve_switch)
        {
          if(opts_init.nx != 0) distmem_real_vctrs.insert({&up, 0});
          if(opts_init.ny != 0) distmem_real_vctrs.insert({&vp, 0});
          if(opts_init.nz != 0) distmem_real_vctrs.insert({&wp, 0});
        }

        if(opts_init.turb_cond_switch)
        {
          distmem_real_vctrs.insert({&wp, 0});
          distmem_real_vctrs.insert({&ssp, 0});
          distmem_real_vctrs.insert({&dot_ssp, 0});
        }
         
        if(opts_init.diag_incloud_time)
          distmem_real_vctrs.insert({&incloud_time, detail::no_initial_value});

        // initializing distmem_n_vctrs - list of n_t vectors with properties of SDs that have to be copied/removed/recycled when a SD is copied/removed/recycled
        distmem_n_vctrs.insert(&n);

        // init number of temporary real vctrs
        if(opts_init.chem_switch || allow_sstp_cond || n_dims >= 2)
          tmp_device_real_part.add_vector();
        if((allow_sstp_cond && opts_init.exact_sstp_cond) || n_dims==3 || opts_init.turb_cond_switch)
          tmp_device_real_part.add_vector();
        if(allow_sstp_cond && opts_init.exact_sstp_cond)
        {
          tmp_device_real_part.add_vector();
          tmp_device_real_part.add_vector();
          if(opts_init.const_p)
            tmp_device_real_part.add_vector();
        }

        resize_size_vctrs.insert(&ijk);
        resize_size_vctrs.insert(&sorted_ijk);
        resize_size_vctrs.insert(&sorted_id);
        resize_size_vctrs.insert(&tmp_device_size_part);
        if (opts_init.nx != 0) resize_size_vctrs.insert(&i);
        if (opts_init.ny != 0) resize_size_vctrs.insert(&j);
        if (opts_init.nz != 0) resize_size_vctrs.insert(&k);
      }

      void sanity_checks();
      void init_SD_with_distros();
      void init_SD_with_distros_sd_conc(const common::unary_function<real_t> &, const real_t &);
      void init_SD_with_distros_tail(const common::unary_function<real_t> &, const real_t);
      void init_SD_with_distros_const_multi(const common::unary_function<real_t> &);
      void init_SD_with_distros_finalize(const real_t &, const bool unravel_ijk = true);
      void init_SD_with_sizes();
      void init_sanity_check(
        const arrinfo_t<real_t>, const arrinfo_t<real_t>, const arrinfo_t<real_t>,
        const arrinfo_t<real_t>, const arrinfo_t<real_t>,
        const arrinfo_t<real_t>, const arrinfo_t<real_t>,
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
      void init_n_dry_sizes(const real_t &conc, const thrust_size_t &sd_conc);

      void dist_analysis_sd_conc(
        const common::unary_function<real_t> &n_of_lnrd,
        const n_t sd_conc,
        const real_t dt = 1.
      );
      void dist_analysis_const_multi(
        const common::unary_function<real_t> &n_of_lnrd 
      );
      void reserve_hskpng_npart();
      void init_ijk();
      void init_xyz();
      void init_kappa(const real_t &);
      void init_incloud_time();
      void init_count_num_sd_conc(const real_t & = 1);
      void init_count_num_const_multi(const common::unary_function<real_t> &);
      void init_count_num_const_multi(const common::unary_function<real_t> &, const thrust_size_t &);
      void init_count_num_dry_sizes(const std::pair<real_t, int> &);
      void init_count_num_hlpr(const real_t &, const thrust_size_t &);
      void init_count_num_src(const thrust_size_t &);
      template <class arr_t>
      void conc_to_number(arr_t &arr); 
      void init_e2l(const arrinfo_t<real_t> &, thrust_device::vector<real_t>*, const int = 0, const int = 0, const int = 0, const long int = 0);
      void init_wet();
      void init_sync();
      void init_grid();
      void init_hskpng_ncell();
      void init_chem();
      void init_chem_aq();
      void init_sstp();
      void init_sstp_chem();
      void init_kernel();
      void init_vterm();

      void fill_outbuf(thrust::host_vector<real_t>&);
      std::vector<real_t> fill_attr_outbuf(const std::string&);
      void mpi_exchange();

           // rename hskpng_ -> step_?
      void hskpng_sort_helper(bool);
      void hskpng_sort();
      void hskpng_shuffle_and_sort();
      void hskpng_count();
      void ravel_ijk(const thrust_size_t begin_shift = 0);
      void unravel_ijk(const thrust_size_t begin_shift = 0);
      void hskpng_ijk();
      void hskpng_Tpr();
      void hskpng_mfp();

      void hskpng_vterm_all();
      void hskpng_vterm_invalid();
      void hskpng_tke();
      void hskpng_turb_vel(const real_t &dt, const bool only_vertical = false);
      void hskpng_turb_dot_ss();
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
        const typename thrust_device::vector<real_t>::iterator &vec_bgn,
        const thrust_size_t npart,
        const bool cons
      ); 
      void moms_rng(
        const real_t &min, const real_t &max, 
        const typename thrust_device::vector<real_t>::iterator &vec_bgn,
        const bool cons
      ); 
      void moms_calc(
        const typename thrust_device::vector<real_t>::iterator &vec_bgn,
        const thrust_size_t npart,
        const real_t power,
        const bool specific = true
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
        arrinfo_t<real_t> &
      );

      void adjust_timesteps(const real_t &dt);
      void adve();
      void turb_adve(const real_t &dt);
      template<class adve_t>
      void adve_calc(bool, thrust_size_t = 0);
      void sedi(const real_t &dt);
      void subs(const real_t &dt);

      void cond_dm3_helper();
      void cond(const real_t &dt, const real_t &RH_max, const bool turb_cond);
      void cond_sstp(const real_t &dt, const real_t &RH_max, const bool turb_cond);
      template<class pres_iter, class RH_iter>
      void cond_sstp_hlpr(const real_t &dt, const real_t &RH_max, const thrust_device::vector<real_t> &Tp, const pres_iter &pi, const RH_iter &rhi);
      void update_th_rv(thrust_device::vector<real_t> &);
      void update_state(thrust_device::vector<real_t> &, thrust_device::vector<real_t> &);
      void update_pstate(thrust_device::vector<real_t> &, thrust_device::vector<real_t> &);
      void update_incloud_time(const real_t &dt);

      void coal(const real_t &dt, const bool &turb_coal);

      void chem_vol_ante();
      void chem_flag_ante();
      void chem_henry(const real_t &dt);
      void chem_dissoc();
      void chem_react(const real_t &dt);
      void chem_cleanup();
      void chem_post_step();
 
      thrust_size_t rcyc();
      void bcnd();

      void src(const real_t &dt);
      void src_dry_distros_simple(const real_t &dt);
      void src_dry_distros_matching(const real_t &dt);
      void src_dry_distros(const real_t &dt);
      void src_dry_sizes(const real_t &dt);

      void rlx(const real_t);
      void rlx_dry_distros(const real_t);

      void sstp_step(const int &step);
      void sstp_step_exact(const int &step);
      void sstp_step_ssp(const real_t &dt);
      void sstp_save();
      void sstp_step_chem(const int &step);
      void sstp_save_chem();

      void post_copy(const opts_t<real_t>&);

      // two functions for calculating changes in rv and th due to condensation on SDs initialized during simulation, e.g. via source or relaxation
      // NOTE: curently not used, because of small sizes of these droplets
      void ante_adding_SD();
      void post_adding_SD();

      // distmem stuff
      void xchng_domains();
      void xchng_courants();
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
