#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string 
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_blk_1m";

  auto n = h5n(h5);

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    for (auto &plt : std::set<std::string>({"rc", "rr"}))
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

      if (plt == "rr")
      {
	auto rr = h5load(h5, "rr", at * n["outfreq"]) * 1e3;
	gp << "set logscale cb\n";
	gp << "set title 'rain water mixing ratio r_r [g/kg]'\n";
	gp << "set cbrange [1e-2:1]\n";
	plot(gp, rr);
      }
    }
  }
}
