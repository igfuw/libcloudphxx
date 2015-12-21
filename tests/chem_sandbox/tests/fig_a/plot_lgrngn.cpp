#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_lgrngn";

  auto n = h5n(h5);

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    for (auto &plt : std::set<std::string>({"rl", "rr", "nc", "nr", "ef", "na"}))
    {
      Gnuplot gp;
      init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".svg", 1, 1, n); 

      if (n["x"] == 76 && n["z"] == 76)
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
      else assert(false);
    } // var loop
  } // time loop
} // main
