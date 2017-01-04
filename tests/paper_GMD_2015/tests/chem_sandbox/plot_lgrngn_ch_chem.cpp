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

  std::string dir = string(av[1]) + "/tests/chem_sandbox/";//1_no_chem_results/";
  std::string h5  = dir + "out_hall_pinsky_stratocumulus";

  auto n = h5n(h5);

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    for (auto &plt : std::set<std::string>({"SO2g",  "O3g",    "H2O2g",  "CO2g",   "NH3g",  "HNO3g", 
                                           "S_IV_aq", "S_VI_aq", "O3_aq",  "H2O2_aq", "H_aq",  
                                           "C_IV_aq", "N_III_aq", "N_V_aq"}))

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

      if (plt == "SO2g") //200e-12
      {                                                      // TODO this is lazy (assumes pd = p), do it better
        auto chem = h5load(h5, plt, at * n["outfreq"]) * (M_d<float>() / M_SO2<float>()) * 1e9;
        gp << "set title 'gas vol conc [ppb]'\n";
//        gp << "set cbrange [0.05:0.2]\n";
        plot(gp, chem);
      }
      else if (plt == "O3g") //50e-9
      {
        auto chem = h5load(h5, plt, at * n["outfreq"]) * (M_d<float>() / M_O3<float>()) * 1e9;
        gp << "set title 'gas vol conc [ppb]'\n";
        gp << "set cbrange [50.46:50.62]\n";
        plot(gp, chem);
      }
      else if (plt == "H2O2g") //500e-12 ; 0
      {
        auto chem = h5load(h5, plt, at * n["outfreq"]) * (M_d<float>() / M_H2O2<float>()) * 1e9;
        gp << "set title 'gas vol conc [ppb]'\n";
        gp << "set cbrange [0:0.7]\n";
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
  //      gp << "set cbrange [0:0.1]\n";
        plot(gp, chem);
      }
      else if (plt == "HNO3g") // 100e-12
      {
        auto chem = h5load(h5, plt, at * n["outfreq"]) * (M_d<float>() / M_HNO3<float>()) * 1e9;
        gp << "set title 'gas vol conc [ppb]'\n";
    //    gp << "set cbrange [0:0.1]\n";
        plot(gp, chem);
      }


      else if (plt == "S_IV_aq")
      {
        auto chem = h5load(h5, "chem_S_IV_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'S_IV_aq [ug/kg]'\n";
        //gp << "set cbrange [0:1e-5]\n";
        plot(gp, chem);
      }
      else if (plt == "S_VI_aq")
      {
        auto chem = h5load(h5, "chem_S_VI_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'S_VI_aq [ug/kg]'\n";
      //  gp << "set cbrange [0:1]\n";
        plot(gp, chem);
      }

      else if (plt == "O3_aq")
      {
        auto chem = h5load(h5, "chem_O3_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'O3_aq [ug/kg]'\n";
        gp << "set cbrange [-1e-5:2.8e-5]\n";
        plot(gp, chem);
      }
       else if (plt == "H2O2_aq")
      {
        auto chem = h5load(h5, "chem_H2O2_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'H2O2_aq [ug/kg]'\n";
        gp << "set cbrange [0:1]\n";
        plot(gp, chem);
      }

      else if (plt == "H_aq")
      {
        auto chem = h5load(h5, "chem_H_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'H_aq [ug/kg]'\n";
  //      gp << "set cbrange [0:0.016]\n";
        plot(gp, chem);
      }

       else if (plt == "C_IV_aq")
      {
        auto chem = h5load(h5, "chem_C_IV_aq", at * n["outfreq"]) * 1e9;
        gp << "set title 'C_IV_aq [ug/kg]'\n";
        //gp << "set cbrange [0:2.5]\n";
        plot(gp, chem);
      }

       else if (plt == "N_III_aq")
      {
        auto chem = h5load(h5, "chem_N_III_aq", at * n["outfreq"]) * 1e9; //TODO
        gp << "set title 'N_III_aq [ug/kg]'\n";
        //gp << "set cbrange [0:7e-5]\n";
        plot(gp, chem);
      }

       else if (plt == "N_V_aq")
      {
        auto chem = h5load(h5, "chem_N_V_aq", at * n["outfreq"]) * 1e9; 
        gp << "set title 'N_V_aq [ug/kg]'\n";
  //      gp << "set cbrange [0:128]\n";
        plot(gp, chem);
      }
 
      else assert(false);
    } // var loop
  } // time loop
} // main

