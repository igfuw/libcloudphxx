#include <cstdlib> // system()
#include <set>
#include <string>
#include <sstream> // std::ostringstream

#include <fstream>
#include <iostream>

#include "../common.hpp"
#include "bins.hpp"

using std::ostringstream;
using std::set;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string bins_wet_str, bins_dry_str, bins_wet_str_2;

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

  {
    ostringstream tmp;
    vector<quantity<si::length>> left_edges = bins_wet();
    for (int i = 0; i < left_edges.size()-1; ++i)
      tmp << float(left_edges[i] / si::metres) << ":" << float(left_edges[i + 1] / si::metres) << "|0,3;";
    bins_wet_str_2 = tmp.str();
  }

  {
    string cmn = 
      //"--outfreq=11800  --dt=1 --nt=11800 --spinup=10000 --nx=76 --nz=76 --relax_th_rv=false --rng_seed=44 ";
      "--outfreq=200  --dt=1 --nt=11800 --spinup=10000 --nx=76 --nz=76 --relax_th_rv=false  "
      "--backend=CUDA --adv_serial=False --sd_conc=256 --sstp_cond=10 --coal=True --sedi=True "
      "--w_max=.6 --terminal_velocity=beard77fast ";

    string stats =         
        " --out_wet=\""
          "1.e-6:25e-6|0,1,2,3;"   // "cloud"
          "25e-6:1|0,3;"           // "rain"
          + bins_wet_str_2 +       // aerosol spectrum (wet)
        "\""
        " --out_dry=\""
          "0.:1.|0,1;"
          + bins_dry_str +         // aerosol spectrum (dry)
        "\"";

    string chem_stats =         
        " --out_chem=\""
          "0:1|0;"                 // chem spectrum (dry)
        "\""
        " --out_wet_pH=\""
          + bins_wet_str +         // spectrum for S_VI and H+ (wet)
        "\"";  

    string size1 = "--mean_rd1=0.05e-6 --sdev_rd1=1.8 --n1_stp=50e6  --mean_rd2=0.1e-6 --sdev_rd2=1.5 --n2_stp=0 ";
    string size2 = "--mean_rd1=0.05e-6 --sdev_rd1=1.8 --n1_stp=150e6 --mean_rd2=0.1e-6 --sdev_rd2=1.5 --n2_stp=0 ";
    string size3 = "--mean_rd1=0.05e-6 --sdev_rd1=1.8 --n1_stp=0     --mean_rd2=0.1e-6 --sdev_rd2=1.5 --n2_stp=50e6 ";
    string size4 = "--mean_rd1=0.05e-6 --sdev_rd1=1.8 --n1_stp=0     --mean_rd2=0.1e-6 --sdev_rd2=1.5 --n2_stp=150e6 ";

    //chem_rct switched on afetr spinup in set_chem 
    string chem_cmn = "--micro=lgrngn_chem --kernel=hall_pinsky_stratocumulus --rng_seed=13 "
                      "--chem_switch=True --chem_dsl=True --chem_dsc=True --chem_rho=1.8e3 --sstp_chem=10 ";

    string ch_ini_base   = "--SO2_g_0=.2e-9 --O3_g_0=25e-9 --H2O2_g_0=.4e-9 --CO2_g_0=360e-6 --NH3_g_0=.1e-9 --HNO3_g_0=.1e-9 ";
    string ch_ini_case3  = "--SO2_g_0=.2e-9 --O3_g_0=25e-9 --H2O2_g_0=.4e-9 --CO2_g_0=360e-6 --NH3_g_0=.4e-9 --HNO3_g_0=.1e-9 ";
    string ch_ini_case1a = "--SO2_g_0=.2e-9 --O3_g_0=0     --H2O2_g_0=0     --CO2_g_0=360e-6 --NH3_g_0=.1e-9 --HNO3_g_0=.1e-9 ";

    string case1   = "--micro=lgrngn --outdir=out_case1 --kernel=hall_pinsky_stratocumulus --rng_seed=13 ";
    string case1a  = "--outdir=out_case1a ";
    string caseb   = "--outdir=out_case_base ";
    string case3   = "--outdir=out_case3 ";
    string case4   = "--outdir=out_case4 ";
    string case5   = "--outdir=out_case5 ";
    string case6   = "--outdir=out_case6_no_coll_no_sedi ";

    set<string> opts_micro({
        //cmn + " " + case1  + " " + stats + " " + size1,                                                           // case1
        //cmn + " " + chem_cmn + " " + stats + " " + chem_stats + " " + size1 + " " + ch_ini_base   + " " + caseb,  // case base
        //cmn + " " + chem_cmn + " " + stats + " " + chem_stats + " " + size1 + " " + ch_ini_case1a + " " + case1a, // case1a
        //cmn + " " + chem_cmn + " " + stats + " " + chem_stats + " " + size1 + " " + ch_ini_case3  + " " + case3,  // case3
        //cmn + " " + chem_cmn + " " + stats + " " + chem_stats + " " + size2 + " " + ch_ini_base   + " " + case4,  // case4
        //cmn + " " + chem_cmn + " " + stats + " " + chem_stats + " " + size3 + " " + ch_ini_base   + " " + case5,  // case5
        cmn + " " + chem_cmn + " " + stats + " " + chem_stats + " " + size4 + " " + ch_ini_base   + " " + case6   // case6
    });

    // run the above simulation cases
    for (auto &opts_m : opts_micro)
    {
      ostringstream cmd;
      cmd << av[1] << "/src/icicle " << opts_m;  
      notice_macro("about to call: " << cmd.str())

      std::ofstream log;
      log.open("command_log", std::ios::app);
      log << cmd.str() << "\n";
      log.close();  

      if (EXIT_SUCCESS != system(cmd.str().c_str()))
        error_macro("model run failed: " << cmd.str())
    }
/*
    // run case2 sumulations (all collision kernels and 5 different random seeds)
    for (const std::string kernel : {"hall", "hall_davis_no_waals", "hall_pinsky_stratocumulus", 
                                     "onishi_hall", "onishi_hall_davis_no_waals", "vohl_davis_no_waals"})
    {
      for (const std::string seed : {"13", "42", "44", "9", "6"})
      {
        string case2 = cmn+ " " + "--micro=lgrngn --rng_seed="+seed+ " --outdir=out_case2_"+kernel+"_"+seed+ " --kernel="+kernel + " " +stats+ " " +size1;

        ostringstream cmd;
        cmd << av[1] << "/src/icicle " << case2;  
        notice_macro("about to call: " << cmd.str())

        std::ofstream log;
        log.open("command_log", std::ios::app);
        log << cmd.str() << "\n";
        log.close();  

        if (EXIT_SUCCESS != system(cmd.str().c_str()))
          error_macro("model run failed: " << cmd.str())
      }
    }*/
  }
}
