#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/moist_air.hpp>

using namespace libcloudphxx::common::molar_mass;
using namespace libcloudphxx::common::moist_air;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string dir = string(av[1]) + "/tests/chem_sandbox/";

  for (const std::string sim_run :
           {
            "out_case1",
            "out_case2_hall",
            "out_case2_hall_davis_no_waals",
            "out_case2_hall_pinsky_stratocumulus",
            "out_case2_onishi_hall",
            "out_case2_onishi_hall_davis_no_waals",
            "out_case2_vohl_davis_no_waals",
            "out_case_base",
            "out_case1a",
            "out_case3",
            "out_case4",
            "out_case5",
            "out_case6"
           })

  {
    std::string h5  = dir + sim_run;
  
    auto n = h5n(h5);
    for (int at = 0; at < n["t"]; ++at)
    {
      if (at * n["outfreq"] >= 10000 || at == 0) //plot only after spinup
      {
        for (auto &plt : std::set<std::string>({"rl", "rr", "nc", "nr", "ef", "na", "rd", "sd_conc", "th", "rv"})) 
        {
          Gnuplot gp;
          init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".svg", 1, 1, n); 

          if (plt == "rl")
          {
  	    // cloud water content
  	    //                                                         rho_w  kg2g
  	    auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
  	    gp << "set title 'cloud water mixing ratio [g/kg]'\n";
  	    gp << "set cbrange [0:1.6]\n";
  	    plot(gp, tmp);
          }
  
          else if (plt == "rr")
          {
  	    // rain water content
  	    //                                                         rho_w  kg2g
  	    auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
  	    gp << "set logscale cb\n";
  	    gp << "set title 'rain water mixing ratio [g/kg]'\n";
  	    gp << "set cbrange [1e-2:1]\n";
  	    plot(gp, tmp);
  	    gp << "unset logscale cb\n";
          }
  
          else if (plt == "nc")
          {
  	    // cloud particle concentration
  	    auto tmp = 1e-6 * h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
  	    gp << "set title 'cloud droplet spec. conc. [mg^{-1}]'\n";
  	    gp << "set cbrange [0:160]\n";
  	    plot(gp, tmp);
          }
  
          else if (plt == "nr")
          {
  	    // rain particle concentration
  	    auto tmp = 1e-6 * h5load(h5, "rw_rng001_mom0", at * n["outfreq"]);
  	    gp << "set title 'rain drop spec. conc. [mg^{-1}]'\n";
  	    gp << "set cbrange [.01:10]\n";
  	    gp << "set logscale cb\n";
  	    plot(gp, tmp);
  	    gp << "unset logscale cb\n";
          }
  
          else if (plt == "ef")
          {
  	    // effective radius
  	    auto r_eff = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) / h5load(h5, "rw_rng000_mom2", at * n["outfreq"]) * 1e6;
  	    gp << "set title 'cloud droplet effective radius [μm]'\n"; 
  	    gp << "set cbrange [1:14]\n";
  	    plot(gp, r_eff);
          }
  
          else if (plt == "na")
          {
      	    // aerosol concentration
    	    blitz::Array<float, 2> tmp(h5load(h5, "rw_rng002_mom0", at * n["outfreq"]));
    	    vector<quantity<si::length>> left_edges = bins_wet();
    	    for (int i = 1; i < left_edges.size()-1; ++i)
     	    {
    	      if (left_edges[i + 1] > 1e-6 * si::metres) break;
    	      ostringstream str;
    	      str << "rw_rng" << std::setw(3) << std::setfill('0') << i + 2  << "_mom0";
    	      tmp = tmp + h5load(h5, str.str(), at * n["outfreq"]);
    	    }
    	    gp << "set cbrange [" << 0 << ":" << 160 << "]\n";
    	    gp << "set title 'aerosol concentration [mg^{-1}]'\n";
    	    tmp /= 1e6;
    	    plot(gp, tmp);
          }
    
          else if (plt == "rd")
          {
            auto r_d = h5load(h5, "rd_rng000_mom1", at * n["outfreq"])/h5load(h5, "rd_rng000_mom0", at * n["outfreq"]) * 1e6;
    	    gp << "set title 'dry radius [μm]'\n"; 
     	    gp << "set cbrange [0:0.14]\n";
     	    //gp << "set cbrange [0:0.1]\n";
    	    plot(gp, r_d);
          }
  
          else if (plt == "sd_conc")
          {
            // super-droplet concentration
            auto sd_conc = h5load(h5, "sd_conc", at * n["outfreq"]);
            gp << "set title 'super-droplet concentration [dv-1]'\n";
            gp << "set cbrange [0:256]\n";
            plot(gp, sd_conc);
          }
    
          else if (plt == "th")
          {
            auto th = h5load(h5, "th", at * n["outfreq"]);
            gp << "set title 'potential temperature [K]'\n";
            gp << "set cbrange [289.5:292.5]\n";
            plot(gp, th);
          }
  
          else if (plt == "rv")
          {
            auto rv = h5load(h5, "rv", at * n["outfreq"]) * 1000;
            gp << "set title 'water vapour mixing ratio [g/kg]'\n";
            gp << "set cbrange [6.5:7.5]\n";
            plot(gp, rv);
          }
  
          else assert(false);
        } // var loop
      } // if for only after spinup
    } // time loop
  } // sim loop
} // main

