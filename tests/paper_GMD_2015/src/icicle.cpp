/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/bcond/open_2d.hpp>
#include <libmpdata++/concurr/boost_thread.hpp> // not to conflict with OpenMP used via Thrust in libcloudph++
#include <libmpdata++/concurr/serial.hpp> // not to conflict with OpenMP used via Thrust in libcloudph++

#include "icmw8_case1.hpp" // 8th ICMW case 1 by Wojciech Grabowski)
namespace setup = icmw8_case1;

#include "opts_blk_1m.hpp"
#include "opts_blk_2m.hpp"
#include "opts_lgrngn.hpp"

#include "panic.hpp"

// model run logic - the same for any microphysics
template <class solver_t>
void run(int nx, int nz, int nt, const std::string &outdir, const int &outfreq, int spinup, bool serial, bool relax_th_rv)
{
  // instantiation of structure containing simulation parameters
  typename solver_t::rt_params_t p;

  // output and simulation parameters
  p.grid_size = {nx, nz};
  p.outdir = outdir;
  p.outfreq = outfreq;
  p.spinup = spinup;
  p.relax_th_rv = relax_th_rv;
  setup::setopts(p, nx, nz);
  setopts_micro<solver_t>(p, nx, nz, nt);

  // solver instantiation
  std::unique_ptr<
    concurr::any<
      typename solver_t::real_t, 
      solver_t::n_dims
    >
  > slv;
  if (serial)
  {
    using concurr_t = concurr::serial<
      solver_t, 
      bcond::cyclic, bcond::cyclic,
      bcond::open,   bcond::open 
    >;
    slv.reset(new concurr_t(p));

    // initial condition
    setup::intcond(*static_cast<concurr_t*>(slv.get()));
  }
  else
  {
    using concurr_t = concurr::boost_thread<
      solver_t, 
      bcond::cyclic, bcond::cyclic,
      bcond::open,   bcond::open 
    >;
    slv.reset(new concurr_t(p));

    // initial condition
    setup::intcond(*static_cast<concurr_t*>(slv.get()));
  }


  // setup panic pointer and the signal handler
  panic = slv->panic_ptr();
  set_sigaction();
 
  // timestepping
  slv->advance(nt);
}


// libmpdata++'s compile-time parameters
struct ct_params_common : ct_params_default_t
{
  using real_t = setup::real_t;
  enum { n_dims = 2 };
  enum { opts = opts::nug | opts::fct }; 
  enum { rhs_scheme = solvers::euler_b };
};


// all starts here with handling general options 
int main(int argc, char** argv)
{
  // making argc and argv global
  ac = argc;
  av = argv;

  {
    // note: all options should have default values here to make "--micro=? --help" work
    opts_main.add_options()
      ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
      ("nx", po::value<int>()->default_value(76) , "grid cell count in horizontal")
      ("nz", po::value<int>()->default_value(76) , "grid cell count in vertical")
      ("nt", po::value<int>()->default_value(3600) , "timestep count")
      ("outdir", po::value<std::string>(), "output file name (netCDF-compatible HDF5)")
      ("outfreq", po::value<int>(), "output rate (timestep interval)")
      ("spinup", po::value<int>()->default_value(2400) , "number of initial timesteps during which rain formation is to be turned off")
      ("adv_serial", po::value<bool>()->default_value(false), "force advection to be computed on single thread")
      ("relax_th_rv", po::value<bool>()->default_value(true) , "relaxation of th and rv")
      ("help", "produce a help message (see also --micro X --help)")
    ;
    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); // ignores unknown

    // hendling the "help" option
    if (ac == 1 || (vm.count("help") && !vm.count("micro"))) 
    {
      std::cout << opts_main;
      exit(EXIT_SUCCESS);
    }

    // checking if all required options present
    po::notify(vm); 
    
    // handling outdir && outfreq
    std::string outdir; 
    int outfreq;
    if (!vm.count("help"))
    {
      if (!vm.count("outdir")) throw po::required_option("outdir");
      if (!vm.count("outfreq")) throw po::required_option("outfreq");
      outdir = vm["outdir"].as<std::string>();
      outfreq = vm["outfreq"].as<int>();
    }

    // handling nx, nz, nt options
    int 
      nx = vm["nx"].as<int>(),
      nz = vm["nz"].as<int>(),
      nt = vm["nt"].as<int>(),
      spinup = vm["spinup"].as<int>();

    // handling serial-advection-forcing flag
    bool adv_serial = vm["adv_serial"].as<bool>();

    // handling relaxation flag
    bool relax_th_rv = vm["relax_th_rv"].as<bool>();

    // handling the "micro" option
    std::string micro = vm["micro"].as<std::string>();

    if (micro == "blk_1m")
    {
      // libmpdata++'s compile-time parameters
      if (relax_th_rv)
      {
        struct ct_params_t : ct_params_common
        {
  	  enum { n_eqns = 4 };
          struct ix { enum {th, rv, rc, rr}; };
        };
        run<kin_cloud_2d_blk_1m<ct_params_t>>(nx, nz, nt, outdir, outfreq, spinup, adv_serial, relax_th_rv);
      }
      else
      {
        struct ct_params_t : ct_params_common
        {
          enum { n_eqns = 4 };
          struct ix { enum {th, rv, rc, rr}; };
          enum { hint_norhs = opts::bit(ix::th) | opts::bit(ix::rv) };
        };
        run<kin_cloud_2d_blk_1m<ct_params_t>>(nx, nz, nt, outdir, outfreq, spinup, adv_serial, relax_th_rv);
      }
    }

    else
    if (micro == "blk_2m")
    {
      struct ct_params_t : ct_params_common
      {
	enum { n_eqns = 6 };
	struct ix { enum {th, rv, rc, rr, nc, nr}; }; 
      };
      run<kin_cloud_2d_blk_2m<ct_params_t>>(nx, nz, nt, outdir, outfreq, spinup, adv_serial, relax_th_rv);
    }

    else 
    if (micro == "lgrngn")
    {
      if (relax_th_rv)
      {
        struct ct_params_t : ct_params_common
        {
  	  enum { n_eqns = 2 };
  	  struct ix { enum {th, rv}; };
        };
        run<kin_cloud_2d_lgrngn<ct_params_t>>(nx, nz, nt, outdir, outfreq, spinup, adv_serial, relax_th_rv);
      }
      else
      {
        struct ct_params_t : ct_params_common
        {
  	  enum { n_eqns = 2 };
  	  struct ix { enum {th, rv}; };
          enum { hint_norhs = opts::bit(ix::th) | opts::bit(ix::rv) };
        };
        run<kin_cloud_2d_lgrngn<ct_params_t>>(nx, nz, nt, outdir, outfreq, spinup, adv_serial, relax_th_rv);
      }
    }
    else throw
      po::validation_error(
        po::validation_error::invalid_option_value, micro, "micro" 
      );
  }
}
