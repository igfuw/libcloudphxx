#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include "../../../../include/libcloudph++/common/molar_mass.hpp"
#include "../../../../include/libcloudph++/common/moist_air.hpp"

using namespace libcloudphxx::common::molar_mass;
using namespace libcloudphxx::common::moist_air;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_lgrngn";

  auto n = h5n(h5);

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    for (auto &plt : std::set<std::string>({"rl", "rr", "nc", "nr", "ef", "na", "rd", "sd_conc", "th", "rv", 
                                           "SO2g",  "O3g",    "H2O2g",  "CO2g",   "NH3g",  "HNO3g", 
                                           "SO2aq", "HSO3aq", "SO3aq",  "HSO4aq", "SO4aq", "S_VIaq", 
                                           "O3aq",  "H2O2aq", "Haq",    "OHaq",  
                                           "CO2aq", "HCO3aq", "CO3aq", 
                                           "NH3aq", "NH4aq",  "HNO3aq", "NO3aq"}))

    {
      Gnuplot gp;
      init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".svg", 1, 1, n); 

      if (n["x"] == 50 && n["z"] == 50) //76
      {
	{
	  char lbl = 'i';
	  for (auto &fcs : std::set<std::set<std::pair<int, int>>>({focus.first, focus.second}))
	  {
	    for (auto &pr : fcs) 
	    {
	      auto &x = pr.first;
	      auto &y = pr.second;

	      // black square
	      gp << "set arrow from " << x-1 << "," << y-1 << " to " << x+2 << "," << y-1 << " nohead lw 4 lc rgbcolor '#ffffff' front\n";
	      gp << "set arrow from " << x-1 << "," << y+2 << " to " << x+2 << "," << y+2 << " nohead lw 4 lc rgbcolor '#ffffff' front\n";
	      gp << "set arrow from " << x-1 << "," << y-1 << " to " << x-1 << "," << y+2 << " nohead lw 4 lc rgbcolor '#ffffff' front\n";
	      gp << "set arrow from " << x+2 << "," << y-1 << " to " << x+2 << "," << y+2 << " nohead lw 4 lc rgbcolor '#ffffff' front\n";
	      // white square
	      gp << "set arrow from " << x-1 << "," << y-1 << " to " << x+2 << "," << y-1 << " nohead lw 2 front\n";
	      gp << "set arrow from " << x-1 << "," << y+2 << " to " << x+2 << "," << y+2 << " nohead lw 2 front\n";
	      gp << "set arrow from " << x-1 << "," << y-1 << " to " << x-1 << "," << y+2 << " nohead lw 2 front\n";
	      gp << "set arrow from " << x+2 << "," << y-1 << " to " << x+2 << "," << y+2 << " nohead lw 2 front\n";

	      lbl -= 2;
	    }
	    lbl = 'j';
	  }
	}

	// labels
	{
	  char lbl = 'i';
	  for (auto &fcs : std::set<std::set<std::pair<int, int>>>({focus.first, focus.second}))
	  {
	    for (auto &pr : fcs) 
	    {
	      auto &x = pr.first;
	      auto &y = pr.second;

	      // labels
	      gp << "set label " << int(lbl) << " '" << lbl << "' at " << x+(((lbl+1)/2)%2?-6:+4) << "," << y+.5 << " front font \",20\"\n";

	      lbl -= 2;
	    }
	    lbl = 'j';
	  }
	}
      }

      if (plt == "rl")
      {
	// cloud water content
	//                                                         rho_w  kg2g
	auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
	gp << "set title 'cloud water mixing ratio [g/kg]'\n";
	gp << "set cbrange [0:1.3]\n";
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
	gp << "set cbrange [0:150]\n";
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
	gp << "set cbrange [1:20]\n";
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
	gp << "set cbrange [" << 0 << ":" << 150 << "]\n";
	gp << "set title 'aerosol concentration [mg^{-1}]'\n";
	tmp /= 1e6;
	plot(gp, tmp);
      }
      else if (plt == "rd")
      {
	auto r_d = h5load(h5, "rd_rng000_mom1", at * n["outfreq"])/h5load(h5, "rd_rng000_mom0", at * n["outfreq"]) * 1e6;
	gp << "set title 'dry radius [μm]'\n"; 
 	gp << "set cbrange [0.03:0.1]\n";
	plot(gp, r_d);
      }

      else if (plt == "sd_conc")
      {
        // super-droplet concentration
        auto sd_conc = h5load(h5, "sd_conc", at * n["outfreq"]);
        gp << "set title 'super-droplet concentration [dv-1]'\n";
        gp << "set cbrange [0:128]\n";
        plot(gp, sd_conc);
      }
      else if (plt == "th")
      {
        auto th = h5load(h5, "th", at * n["outfreq"]);
        gp << "set title 'potential temperature [K]'\n";
        gp << "set cbrange [289.5:292]\n";
        plot(gp, th);
      }
      else if (plt == "rv")
      {
        auto rv = h5load(h5, "rv", at * n["outfreq"]) * 1000;
        gp << "set title 'water vapour mixing ratio [g/kg]'\n";
        gp << "set cbrange [6.5:7.7]\n";
        plot(gp, rv);
      }
      else if (plt == "SO2g") //200e-12
      {                                                      // TODO this is lazy (assumes pd = p), do it better
        auto chem = h5load(h5, plt, at * n["outfreq"]) * (M_d<float>() / M_SO2<float>()) * 1e9;
        gp << "set title 'gas vol conc [ppb]'\n";
        gp << "set cbrange [0.05:0.2]\n";
        plot(gp, chem);
      }
      else if (plt == "O3g") //50e-9
      {
        auto chem = h5load(h5, plt, at * n["outfreq"]) * (M_d<float>() / M_O3<float>()) * 1e9;
        gp << "set title 'gas vol conc [ppb]'\n";
        gp << "set cbrange [50.56:50.6]\n";
        plot(gp, chem);
      }
      else if (plt == "H2O2g") //500e-12 ; 0
      {
        auto chem = h5load(h5, plt, at * n["outfreq"]) * (M_d<float>() / M_H2O2<float>()) * 1e9;
        gp << "set title 'gas vol conc [ppb]'\n";
        gp << "set cbrange [0:0.5]\n";
        plot(gp, chem);
      }
      else if (plt == "CO2g") //360e-6
      {
        auto chem = h5load(h5, plt, at * n["outfreq"]) * (M_d<float>() / M_CO2<float>()) * 1e6;
        gp << "set title 'gas vol conc [ppm]'\n";
        //gp << "set cbrange [360:361.2]\n";
        plot(gp, chem);
      }
      else if (plt == "NH3g") // 100e-12
      {
        auto chem = h5load(h5, plt, at * n["outfreq"]) * (M_d<float>() / M_NH3<float>()) * 1e9;
        gp << "set title 'gas vol conc [ppb]'\n";
        gp << "set cbrange [0:0.1]\n";
        plot(gp, chem);
      }
      else if (plt == "HNO3g") // 100e-12
      {
        auto chem = h5load(h5, plt, at * n["outfreq"]) * (M_d<float>() / M_HNO3<float>()) * 1e9;
        gp << "set title 'gas vol conc [ppb]'\n";
        gp << "set cbrange [0:0.1]\n";
        plot(gp, chem);
      }


      else if (plt == "SO2aq")
      {
        auto sd_conc = h5load(h5, "chem_SO2_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'SO2_aq [ug/kg]'\n";
        gp << "set cbrange [0:1e-5]\n";
        plot(gp, sd_conc);
      }
      else if (plt == "HSO3aq")
      {
        auto sd_conc = h5load(h5, "chem_HSO3_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'HSO3_aq [ug/kg]'\n";
        gp << "set cbrange [0:0.06]\n";
        plot(gp, sd_conc);
      }
      else if (plt == "SO3aq")
      {
        auto sd_conc = h5load(h5, "chem_SO3_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'SO3_aq [ug/kg]'\n";
        gp << "set cbrange [0:0.006]\n";
        plot(gp, sd_conc);
      }
      else if (plt == "HSO4aq")
      {
        auto sd_conc = h5load(h5, "chem_HSO4_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'HSO3_aq [ug/kg]'\n";
        gp << "set cbrange [0:0.7]\n";
        plot(gp, sd_conc);
      }
      else if (plt == "SO4aq")
      {
        auto sd_conc = h5load(h5, "chem_SO4_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'SO4_aq [ug/kg]'\n";
        gp << "set cbrange [0:1]\n";
        plot(gp, sd_conc);
      }
      else if (plt == "S_VIaq")
      {
        auto sd_conc = h5load(h5, "chem_S_VI_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'S_VI_aq [ug/kg]'\n";
        gp << "set cbrange [0:1]\n";
        plot(gp, sd_conc);
      }

      else if (plt == "O3aq")
      {
        auto sd_conc = h5load(h5, "chem_O3_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'O3_aq [ug/kg]'\n";
        gp << "set cbrange [0:2.5e-5]\n";
        plot(gp, sd_conc);
      }
       else if (plt == "H2O2aq")
      {
        auto sd_conc = h5load(h5, "chem_H2O2_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'H2O2_aq [ug/kg]'\n";
        gp << "set cbrange [0:1.2]\n";
        plot(gp, sd_conc);
      }

      else if (plt == "Haq")
      {
        auto sd_conc = h5load(h5, "chem_H_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'H_aq [ug/kg]'\n";
        gp << "set cbrange [0:0.016]\n";
        plot(gp, sd_conc);
      }
       else if (plt == "OHaq")
      {
        auto sd_conc = h5load(h5, "chem_OH_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'OH_aq [ug/kg]'\n";
        gp << "set cbrange [0:1.2e-4]\n";
        plot(gp, sd_conc);
      }

       else if (plt == "CO2aq")
      {
        auto sd_conc = h5load(h5, "chem_CO2_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'CO2_aq [ug/kg]'\n";
        gp << "set cbrange [0:2.5]\n";
        plot(gp, sd_conc);
      }
       else if (plt == "HCO3aq")
      {
        auto sd_conc = h5load(h5, "chem_HCO3_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'HCO3_aq [ug/kg]'\n";
        gp << "set cbrange [0:0.3]\n";
        plot(gp, sd_conc);
      }
       else if (plt == "CO3aq")
      {
        auto sd_conc = h5load(h5, "chem_CO3_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'CO3_aq [ug/kg]'\n";
        gp << "set cbrange [0:6e-6]\n";
        plot(gp, sd_conc);
      }

       else if (plt == "NH3aq")
      {
        auto sd_conc = h5load(h5, "chem_NH3_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'NH3_aq [ug/kg]'\n";
        gp << "set cbrange [0:7e-5]\n";
        plot(gp, sd_conc);
      }
       else if (plt == "NH4aq")
      {
        auto sd_conc = h5load(h5, "chem_NH4_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'NH4_aq [ug/kg]'\n";
        gp << "set cbrange [0:0.22]\n";
        plot(gp, sd_conc);
      }

       else if (plt == "HNO3aq")
      {
        auto sd_conc = h5load(h5, "chem_HNO3_aq", at * n["outfreq"]) * 1e9; 
        gp << "set title 'HNO3_aq [ug/kg]'\n";
  //      gp << "set cbrange [0:128]\n";
        plot(gp, sd_conc);
      }
        else if (plt == "NO3aq")
      {
        auto sd_conc = h5load(h5, "chem_NO3_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'NO3_aq [ug/kg]'\n";
  //      gp << "set cbrange [0:128]\n";
        plot(gp, sd_conc);
      }
 
      else assert(false);
    } // var loop
  } // time loop
} // main
