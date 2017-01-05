#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

#include <map>

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/tests/chem_sandbox/",
    h5  = dir + "out_hall_pinsky_stratocumulus",
    svg = dir + "out_hall_pinsky_stratocumulus_chem.svg";

  Gnuplot gp;

  int off = 1; // TODO!!!
  float ymin = 1e-11, ymax = 1e-2;
  const int at = 9000;//11800;

  gp << "set term svg dynamic enhanced fsize 15 size 900, 1500 \n";
  gp << "set output '" << svg << "'\n";
  gp << "set logscale xy\n";
  gp << "set xrange [.002:100]\n";
  gp << "set yrange [" << ymin << ":" << ymax << "]\n";
  gp << "set ylabel '[mg^{-1} μm^{-1}]'\n"; // TODO: add textual description (PDF?)
  gp << "set grid\n";
  gp << "set nokey\n";

  gp << "set xlabel offset 0,1.5 'particle radius [μm]'\n";
  gp << "set key samplen 1.2\n";
  gp << "set xtics rotate by 65 right (.01, .1, 1, 10, 100) \n";

// TODO: use dashed lines to allow printing in black and white... same in image plots

  assert(focus.first.size() == focus.second.size());
  gp << "set multiplot layout " << focus.first.size() << ",2 columnsfirst upwards\n";

  // focus to the gridbox from where the size distribution is plotted
  char lbl = 'i';
  for (auto &fcs : std::set<std::set<std::pair<int,int>>>({focus.first, focus.second}))
  {
    for (auto it = fcs.begin(); it != fcs.end(); ++it)
    {
      const int &x = it->first, &y = it->second;

      gp << "set label 1 '(" << lbl << ")' at graph -.15, 1.02 font ',20'\n";
      //gp << "set title 'x=" << x << " y=" << y << "'\n";

      std::map<float, float> focus_S6;

      //info on the number and location of histogram edges
      vector<quantity<si::length>> left_edges_rw = bins_wet();
      int nsw = left_edges_rw.size() - 1;

      for (int i = 0; i < nsw - off; ++i)
      {
	const string name = "chem_S_VI_rw_rng" + zeropad(i + off) + "_mom0";
	blitz::Array<float, 2> tmp_w(1e6 * h5load(h5, name, at));

	focus_S6[left_edges_rw[i] / 1e-6 / si::metres] = sum(tmp_w(
	  blitz::Range(x-1, x+1),
	  blitz::Range(y-1, y+1)
	)) 
	/ 9 
	/ ((left_edges_rw[i+1] - left_edges_rw[i]) / 1e-6 / si::metres); // per micrometre
      }

      notice_macro("setting-up plot parameters");
      gp << "plot"
	 << "'-' with histeps title 'S6 mass' lw 3 lc rgb 'red' " << endl;
      gp.send(focus_S6);

      lbl -= 2;
    }
    lbl = 'j';
  }
}
