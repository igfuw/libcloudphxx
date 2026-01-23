#include <libcloudph++/lgrngn/factory.hpp>
#include <boost/assign/ptr_map_inserter.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>
#include <iostream>
#include "mpi.h"
#include <numeric>


using namespace std;
using namespace libcloudphxx::lgrngn;
  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;

  //aerosol bimodal lognormal dist. 
  const quantity<si::length, double>
    mean_rd1 = double(30e-6) * si::metres,
    mean_rd2 = double(40e-6) * si::metres;
  const quantity<si::dimensionless, double>
    sdev_rd1 = double(1.4),
    sdev_rd2 = double(1.6);
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, double>
    n1_stp = double(60e6) / si::cubic_metres,
    n2_stp = double(40e6) / si::cubic_metres;


// lognormal aerosol distribution
template <typename T>
struct log_dry_radii : public libcloudphxx::common::unary_function<T>
{
  T funval(const T lnrd) const
  {   
    return T(( 
        lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, double>(lnrd)) +
        lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, double>(lnrd)) 
      ) * si::cubic_metres
    );  
  }   
};  


void two_step(particles_proto_t<double> *prtcls, 
             arrinfo_t<double> th,
             arrinfo_t<double> rhod,
             arrinfo_t<double> rv,
             arrinfo_t<double> Cx,
             arrinfo_t<double> Cy,
             arrinfo_t<double> Cz  ,
             opts_t<double> opts)
{
    prtcls->step_sync(opts,th,rv,rhod,Cx, Cy, Cz);
    prtcls->step_async(opts);
}

int m1(int a)
{
  return a > 1? a : 1;
}


const int nx_min = 2;

void test(backend_t backend, std::string back_name, int ndims, bool dir, int n_devices) // n_devices - number of GPUs used per node, each has to be controlled by a single process
{
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size = -1;
  MPI_Comm_size( MPI_COMM_WORLD, &size ); 
  if(rank==0)
  {
    std::cerr << "ndims: " << ndims <<  " direction: " << dir << " backend: " << back_name << " n_devices: " << n_devices << std::endl;
  } 
  MPI_Barrier(MPI_COMM_WORLD);

  opts_init_t<double> opts_init;
  opts_init.dt=3.;
  opts_init.sstp_coal = 1; 
  opts_init.kernel = kernel_t::geometric;
  opts_init.terminal_velocity = vt_t::beard76;
  opts_init.dx = 1;
  opts_init.nx = rank+nx_min;// nx_factor/2*(rank/2+1); // previously, two GPUs on the same node had the same nx - why? 
  int nx_total = (nx_min + (nx_min + size - 1)) / 2. * size;
  //opts_init.nx = nx_min;
  //int nx_total = nx_min * size;
  opts_init.x1 = opts_init.nx * opts_init.dx;// nx_factor/2*(rank/2+1);
  opts_init.sd_conc = 64;
  opts_init.n_sd_max = 1000*opts_init.sd_conc;
  opts_init.rng_seed = 4444 + rank;
  if(ndims>1)
  {
    opts_init.dz = 1; 
    opts_init.nz = 3; 
    opts_init.z1 = opts_init.nz * opts_init.dz;
  }
  if(ndims==3)
  {
    opts_init.dy = 1; 
    opts_init.ny = 4; 
    opts_init.y1 = opts_init.ny * opts_init.dy;
  }
  opts_init.dev_id = rank%n_devices; 
  //opts_init.dev_id = rank; 
  std::cout << "device id: " << opts_init.dev_id << std::endl;
//  opts_init.sd_const_multi = 1;

/*
  boost::assign::ptr_map_insert<
    log_dry_radii<double> // value type
  >(  
    opts_init.dry_distros // map
  )(  
    0.001 // key
  ); 
*/

  opts_init.dry_distros.emplace( 
    libcloudphxx::lgrngn::kappa_rd_insol_t<double>{double(0.001), double(0.)}, // kappa, rd_insol
    std::make_shared<log_dry_radii<double>>() // distribution
  );

  particles_proto_t<double> *prtcls;

  prtcls = factory<double>(
    backend,
    opts_init
  );

  std::vector<double> vth(opts_init.nx * m1(opts_init.ny) * opts_init.nz, 300.);
  std::vector<double> vrhod(opts_init.nx * m1(opts_init.ny) * opts_init.nz, 1.225);
  std::vector<double> vrv(opts_init.nx * m1(opts_init.ny) * opts_init.nz, 0.01);
  std::vector<double> vCxm((opts_init.nx + 1) * m1(opts_init.ny) * opts_init.nz, -1);
  std::vector<double> vCxp((opts_init.nx + 1) * m1(opts_init.ny) * opts_init.nz, 1);
  std::vector<double> vCy((opts_init.nx) * (m1(opts_init.ny+1)) * opts_init.nz, 0);
  std::vector<double> vCz((opts_init.nx) * m1(opts_init.ny) * (opts_init.nz+1), 0);

  long int strides[] = {0, 1, 1};
  long int xstrides[] = {0, 1, 1};
  long int ystrides[] = {0, 1, 1};
  long int zstrides[] = {0, 1, 1};

  arrinfo_t<double> th(vth.data(), strides);
  arrinfo_t<double> rhod(vrhod.data(), strides);
  arrinfo_t<double> rv(vrv.data(), strides);
  arrinfo_t<double> Cx( dir ? vCxm.data() : vCxp.data(), xstrides);
  arrinfo_t<double> Cz(vCz.data(), ystrides);
  arrinfo_t<double> Cy(vCy.data(), ystrides);

  if(ndims==1)
    prtcls->init(th,rv,rhod, arrinfo_t<double>(), Cx);
  else if(ndims==2)
    prtcls->init(th,rv, rhod, arrinfo_t<double>(), Cx, arrinfo_t<double>(), Cz);
  else if(ndims==3)
    prtcls->init(th,rv, rhod, arrinfo_t<double>(), Cx, Cy, Cz);

  opts_t<double> opts;
  opts.adve = 0;
  opts.sedi = 0;
  opts.cond = 0;
  opts.coal = 1;
//  opts.chem = 0;

  double *out;

  int n_cell = opts_init.nx * m1(opts_init.ny) * opts_init.nz;
  int n_cell_tot = nx_total * m1(opts_init.ny) * opts_init.nz;
  
  std::vector<int> recvcount(size), displs(size);
  std::iota(recvcount.begin(), recvcount.end(), nx_min); 
  std::transform(recvcount.begin(), recvcount.end(), recvcount.begin(), [opts_init](const double& c){return c*opts_init.nz * m1(opts_init.ny);});
  std::partial_sum(recvcount.begin(), recvcount.end()-1, displs.begin()+1);
  displs[0] = 0;

  MPI_Barrier(MPI_COMM_WORLD);

  // output variables
  std::map<std::string, std::vector<double>> global_out_post_coal,
                                             global_out_post_adve;

  std::vector<std::string> out_names = {"sd_conc", "rd", "rw", "kpa"};

  for(std::string name: out_names)
  {
    global_out_post_coal[name] = std::vector<double>(nx_total * opts_init.nz * m1(opts_init.ny));
    global_out_post_adve[name] = std::vector<double>(nx_total * opts_init.nz * m1(opts_init.ny));
  }
  
  // run the simulation 

  // coalescence time steps
  for(int i=0;i<70;++i)
    two_step(prtcls,th,rhod,rv,Cx, ndims==2 ? arrinfo_t<double>() : Cy, Cz, opts);

  // diagnostics
  for(std::string name: out_names)
  {
    prtcls->diag_all();
    if(name == "sd_conc")
      prtcls->diag_sd_conc();
    else if(name == "rd")
      prtcls->diag_dry_mom(1);
    else if(name == "rw")
      prtcls->diag_wet_mom(1);
    else if(name == "kpa")
      prtcls->diag_kappa_mom(1);
    out = prtcls->outbuf();
    MPI_Gatherv(out, opts_init.nx * opts_init.nz * m1(opts_init.ny), MPI_DOUBLE, global_out_post_coal[name].data(), recvcount.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  // advection time steps
  opts.coal = 0;
  opts.adve = 1;

  for(int i=0; i<nx_total; ++i)
    two_step(prtcls,th,rhod,rv,Cx, ndims==2 ? arrinfo_t<double>() : Cy, Cz, opts);

  // diagnostics
  for(std::string name: out_names)
  {
    prtcls->diag_all();
    if(name == "sd_conc")
      prtcls->diag_sd_conc();
    else if(name == "rd")
      prtcls->diag_dry_mom(1);
    else if(name == "rw")
      prtcls->diag_wet_mom(1);
    else if(name == "kpa")
      prtcls->diag_kappa_mom(1);
    out = prtcls->outbuf();
    MPI_Gatherv(out, opts_init.nx * opts_init.nz * m1(opts_init.ny), MPI_DOUBLE, global_out_post_adve[name].data(), recvcount.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  // test if all variables are the same after advection
  for(std::string name: out_names)
    if(rank==0)
    {
      std::cerr << "global_out_post_coal[" << name <<"]: { ";
      for(auto val: global_out_post_coal[name])
        std::cerr << val << ", ";
      std::cerr << "}" << std::endl;

      std::cerr << "global_out_post_adve[" << name << "]: { ";
      for(auto val: global_out_post_adve[name])
        std::cerr << val << ", ";
      std::cerr << "}" << std::endl;

      if(global_out_post_coal[name] != global_out_post_adve[name])
        throw std::runtime_error("error in advection\n");
    }
}

int main(int argc, char *argv[]){

  int provided_thread_lvl;
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided_thread_lvl);

  std::cout << "hostname: ";
  system("/bin/hostname");

  printf("provided thread lvl: %d\n", provided_thread_lvl);

// parsing arguments

  int n_devices = 1, cuda = 0, opt;
//  int nsecs, tfnd;
//
//  nsecs = 0;
//  tfnd = 0;
//  flags = 0;
  while ((opt = getopt(argc, argv, "c:d:")) != -1) {
      switch (opt) {
      case 'd':
          printf("optarg = %s\n", optarg);
          n_devices = atoi(optarg);
          if(n_devices < 1) throw std::runtime_error("Number of devices (-d option) needs to be greater than 0");
          break;
      case 'c':
          printf("optarg = %s\n", optarg);
          cuda = atoi(optarg);
          break;
      default: /* '?' */
          fprintf(stderr, "Usage: %s [-d number_of_devices_per_node] [-c should cuda be used (bool)]\n",
                  argv[0]);
          exit(EXIT_FAILURE);
      }
  }



  std::vector<backend_t> backends = {backend_t(serial), backend_t(OpenMP)};
  if(cuda)
    backends.push_back(backend_t(CUDA));

  std::map<backend_t, std::string> back_names = {
    { backend_t(serial), "serial" },
    { backend_t(OpenMP), "OpenMP" },
    { backend_t(CUDA), "CUDA" }
  };

  for(auto back: backends)
  {
    // 1d doesnt work with MPI
    // 2D
  MPI_Barrier(MPI_COMM_WORLD);
    test(back, back_names[back], 2, false, n_devices);
  MPI_Barrier(MPI_COMM_WORLD);
    test(back, back_names[back], 2, true, n_devices);
  MPI_Barrier(MPI_COMM_WORLD);
    // 3D
    test(back, back_names[back], 3, false, n_devices);
  MPI_Barrier(MPI_COMM_WORLD);
    test(back, back_names[back], 3, true, n_devices);
  MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Finalize();
}
