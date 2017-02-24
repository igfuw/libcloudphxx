#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_blk_2m";

  auto n = h5n(h5);

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    for (auto &plt : std::set<std::string>({"rc", "rr", "nc", "nr", "th"}))
    {
      Gnuplot gp;
      init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".svg", 1, 1, n); 

      if (plt == "rc") 
      {
	auto rc = h5load(h5, "rc", at * n["outfreq"]) * 1e3;
	gp << "set title 'cloud water mixing ratio r_c [g/kg]'\n"; 
	gp << "set cbrange [0:1.3]\n";
	plot(gp, rc);
      }
      else if (plt == "rr")
      {
	auto rr = h5load(h5, "rr", at * n["outfreq"]) * 1e3;
	gp << "set logscale cb\n";
	gp << "set title 'rain water mixing ratio r_r [g/kg]'\n"; 
	gp << "set cbrange [1e-2:1]\n";
	plot(gp, rr);
	gp << "unset logscale cb\n";
      }
      else if (plt == "nc")
      {
	auto nc = h5load(h5, "nc", at * n["outfreq"]) * 1e-6;
	gp << "set title 'cloud droplet specific concentration n_c [mg^{-1}]'\n"; 
	gp << "set cbrange [0:150]\n";
	plot(gp, nc);
      }
      else if (plt == "nr")
      {
	auto nr = h5load(h5, "nr", at * n["outfreq"]) * 1e-6;
	gp << "set title 'rain drop specific concentration n_r [mg^{-1}]'\n"; 
	gp << "set cbrange [0.01:10]\n";
	gp << "set logscale cb\n";
	plot(gp, nr);
      }
      else if (plt == "th")
      {
	auto nr = h5load(h5, "th", at * n["outfreq"]);
	gp << "set title 'theta_d [K]'\n"; 
	gp << "set cbrange [*:*]\n";
	gp << "unset logscale cb\n";
	plot(gp, nr);
      }
    }
  } 
}
