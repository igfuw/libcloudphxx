#include <cstdlib> // system()
#include <set>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"
#include "bins.hpp"

using std::ostringstream;
using std::set;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string bins_wet_str, bins_dry_str, bins_dry_chem_str;

  {
    ostringstream tmp;
    vector<quantity<si::length>> left_edges = bins_dry();
    for (int i = 0; i < left_edges.size()-1; ++i)
      tmp << float(left_edges[i] / si::metres) << ":" << float(left_edges[i + 1] / si::metres) << "|0;";
    bins_dry_str = tmp.str();
  }

  {
    ostringstream tmp;
    vector<quantity<si::length>> left_edges = bins_wet();
    for (int i = 0; i < left_edges.size()-1; ++i)
      tmp << float(left_edges[i] / si::metres) << ":" << float(left_edges[i + 1] / si::metres) << "|0;";
    bins_wet_str = tmp.str();
  }

  string opts_common = 
  //  "--outfreq=1 --nt=100 --spinup=20 --nx=20 --nz=20 --relax_th_rv=false";
  //    "--outfreq=2000 --nt=90000 --spinup=10   --nx=76 --nz=76 --relax_th_rv=false";
    "--outfreq=200  --nt=8400 --spinup=6600  --nx=20 --nz=20 --relax_th_rv=false";
  set<string> opts_micro({
    "--micro=lgrngn --outdir=out_lgrngn --backend=serial --adv_serial=True --sd_conc=64 --sstp_cond=1 --sstp_chem=10 --coal=False --sedi=False --chem_switch=True --chem_dsl=True --chem_sys_cls=True --chem_dsc=True --chem_rho=1.8e3"  //chem_rct switched on afetr spinup in set_chem  
      " --out_wet=\""
        ".5e-6:25e-6|0,1,2,3;" // FSSP
        "25e-6:1|0,3;"         // "rain"
        + bins_wet_str + // aerosol spectrum (wet)
      "\""
      " --out_dry=\""
        "0.:1.|0,1;"
        + bins_dry_str + // aerosol spectrum (dry)
      "\""
      " --out_chem=\""
        "0:1|0;" // chem spectrum (dry)
      "\""
  });

  for (auto &opts_m : opts_micro)
  {
    ostringstream cmd;
    cmd << av[1] << "/src/icicle " << opts_common << " " << opts_m;  
    notice_macro("about to call: " << cmd.str())

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
  }
}
