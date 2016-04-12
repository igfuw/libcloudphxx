// vim:filetype=cpp
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

#if defined(USE_MPI)
  #include "detail/get_mpi_type.hpp"
#endif

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
      const int n_cell; 
      thrust_size_t n_part,            // total number of SDs
                    n_part_old,        // total number of SDs before source
                    n_part_to_init;    // number of SDs to be initialized by source
      detail::rng<real_t, device> rng;

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
	z;   // z spatial coordinate (for 2D and 3D)

      // dry radii distribution characteristics
      real_t log_rd_min, // logarithm of the lower bound of the distr
             log_rd_max, // logarithm of the upper bound of the distr
             multiplier; // multiplier calculated for the above values

      // terminal velocity (per particle)
      thrust_device::vector<real_t> vt; 
      // sea level term velocity according to Beard 1977, compute once
      thrust_device::vector<real_t> vt_0; 
      // no of bins for cached velocity
      const int vt0_n_bin;
      // ln of min and max radius of cached velocity
      const real_t vt0_ln_r_min, vt0_ln_r_max;

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
        sstp_tmp_rv, // either rv_old or advection-caused change in water vapour mixing ratio
        sstp_tmp_th, // ditto for theta_d
        sstp_tmp_rh, // ditto for rho
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
  
      thrust_device::vector<real_t> 
        T,  // temperature [K]
        p,  // pressure [Pa]
        RH, // relative humisity (p_v / p_vs)
        eta;// dynamic viscosity 

      // sorting needed only for diagnostics and coalescence
      bool sorted;

      // timestep counter
      n_t stp_ctr;

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
      thrust_device::vector<real_t> chem_rhs, chem_ante_rhs, chem_post_rhs;
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
        tmp_device_real_part_chem,  // only allocated if chem_switch==1
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
      unsigned int n_x_bfr;

      // number of cells in devices to the left of this one
      unsigned int n_cell_bfr;

      // x0 of the process to the right
      real_t rgt_x0;

      // x1 of the process to the left
      real_t lft_x1;

      // in/out buffers for SDs copied from other GPUs
      thrust_device::vector<n_t> in_n_bfr, out_n_bfr;
      thrust_device::vector<real_t> in_real_bfr, out_real_bfr;

      // ids of sds to be copied with distmem
      thrust_device::vector<thrust_size_t> &lft_id, &rgt_id;

      // real_t vectors copied in distributed memory case
      std::vector<thrust_device::vector<real_t>*> distmem_real_vctrs;

      // number of real_t vectors to be copied
      const int distmem_real_vctrs_count;


      // methods

      // fills u01[0:n] with random numbers
      void rand_u01(thrust_size_t n) { rng.generate_n(u01, n); }

      // fills un[0:n] with random numbers
      void rand_un(thrust_size_t n) { rng.generate_n(un, n); }

      // max(1, n)
      int m1(int n) { return n == 0 ? 1 : n; }

      // ctor 
      impl(const opts_init_t<real_t> &_opts_init, const std::pair<detail::bcond_t, detail::bcond_t> &bcond, const int &mpi_rank, const int &mpi_size) : 
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
        distmem_real_vctrs_count(
          n_dims == 3 ? 7 :
            n_dims == 2 ? 6 : 
              n_dims == 1 ? 5:
                0),  // distmem doesnt work for 0D anyway
        distmem_real_vctrs(7),
        lft_x1(-1),  // default to no
        rgt_x0(-1),  // MPI boudanry
        lft_id(i),   // note: reuses i vector
        rgt_id(tmp_device_size_part),
        vt0_n_bin(10000),
        vt0_ln_r_min(log(5e-7)),
        vt0_ln_r_max(log(3e-3))  // Beard 1977 is defined on 1um - 6mm diameter range
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
          // check temporarily disabled since dewv_id is not passed anymore, TODO: fix it
//	  if (!(opts_init.x1 > opts_init.x0 && opts_init.x1 <= m1(opts_init.nx) * opts_init.dx) && dev_id == -1) // only for single device runs, since on multi_CUDA x1 is not yet adjusted to local domain
//            throw std::runtime_error("!(x1 > x0 & x1 <= min(1,nx)*dx)");
	  if (!(opts_init.y1 > opts_init.y0 && opts_init.y1 <= m1(opts_init.ny) * opts_init.dy))
            throw std::runtime_error("!(y1 > y0 & y1 <= min(1,ny)*dy)");
	  if (!(opts_init.z1 > opts_init.z0 && opts_init.z1 <= m1(opts_init.nz) * opts_init.dz))
            throw std::runtime_error("!(z1 > z0 & z1 <= min(1,nz)*dz)");
        }

        if (opts_init.dt == 0) throw std::runtime_error("please specify opts_init.dt");
        if (opts_init.sd_conc == 0) throw std::runtime_error("please specify opts_init.sd_conc");
        if (opts_init.coal_switch)
        {
          if(opts_init.terminal_velocity == vt_t::undefined) throw std::runtime_error("please specify opts_init.terminal_velocity or turn off opts_init.coal_switch");
          if(opts_init.kernel == kernel_t::undefined) throw std::runtime_error("please specify opts_init.kernel");
        }
        if (opts_init.sedi_switch)
          if(opts_init.terminal_velocity == vt_t::undefined) throw std::runtime_error("please specify opts_init.terminal_velocity or turn off opts_init.sedi_switch");

        // set 0 dev_count to mark that its not a multi_CUDA spawn
        // if its a spawn, multi_CUDA ctor will alter it
        opts_init.dev_count = 0; 

        // initialising host temporary arrays
        {
          int n_grid;
          switch (n_dims) // TODO: document that 3D is xyz, 2D is xz, 1D is x
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
            case 1:
              n_grid = opts_init.nx+1;
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
        ptr_t arr[] = {&rd3, &rw2, &kpa, &vt, &x, &z, &y};
        distmem_real_vctrs = std::vector<ptr_t>(arr, arr + sizeof(arr) / sizeof(ptr_t) );
      }

      void sanity_checks();

      void init_dry();
      void init_n(
        const real_t kappa, // TODO: map
        const common::unary_function<real_t> *n_of_lnrd
      );
      void dist_analysis(
        const common::unary_function<real_t> *n_of_lnrd,
        const n_t sd_conc,
        const real_t dt = 1.
      );
      void init_ijk();
      void init_xyz();
      void init_count_num();
      void init_e2l(const arrinfo_t<real_t> &, thrust_device::vector<real_t>*, const int = 0, const int = 0, const int = 0, const int = 0);
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
        const real_t power
      );
      void moms_calc_cond(
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
      void update_th_rv(thrust_device::vector<real_t> &);

      void coal(const real_t &dt);

      void chem_vol_ante();
      void chem_flag_ante();
      void chem_henry(const real_t &dt);
      void chem_dissoc();
      void chem_react(const real_t &dt);
      void chem_cleanup();
 
      thrust_size_t rcyc();
      real_t bcnd(); // returns accumulated rainfall

      void src(const real_t &dt);

      void sstp_step(const int &step, const bool &var_rho);
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

    // ctor
    template <typename real_t, backend_t device>
    particles_t<real_t, device>::particles_t(opts_init_t<real_t> opts_init)
    {
      int rank, size;

      // handle MPI init
#if defined(USE_MPI)
      detail::mpi_init(MPI_THREAD_SINGLE, rank, size); 
#else
      rank = 0;
      size = 1;
      // throw an error if ran with mpi, but not compiled for mpi
      if (
        // mpich
        std::getenv("PMI_RANK") != NULL ||
        // openmpi
        std::getenv("OMPI_COMM_WORLD_RANK") != NULL ||
        // lam
        std::getenv("LAMRANK") != NULL ||
        // mvapich2
        std::getenv("MV2_COMM_WORLD_RANK") != NULL
      ) throw std::runtime_error("mpirun environment variable detected but libcloudphxx was compiled with MPI disabled");
#endif
      std::pair<detail::bcond_t, detail::bcond_t> bcond;
      if(size > 1)
        bcond = std::make_pair(detail::distmem_mpi, detail::distmem_mpi);
      else
        bcond = std::make_pair(detail::sharedmem, detail::sharedmem);

      // create impl instance
      pimpl.reset(new impl(opts_init, bcond, rank, size));
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
