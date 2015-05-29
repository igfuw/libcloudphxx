#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

#include <map>

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_lgrngn",
    svg = dir + "out_lgrngn_spec.svg";

  Gnuplot gp;

  int off = 2; // TODO!!!
  float ymin = .4 * .01, ymax = .9 * 10000;
  const int at = 9000;

  gp << "set term svg dynamic enhanced fsize 15 size 900, 1500 \n";
  gp << "set output '" << svg << "'\n";
  gp << "set logscale xy\n";
  gp << "set xrange [.002:100]\n";
  gp << "set yrange [" << ymin << ":" << ymax << "]\n";
  gp << "set ylabel '[mg^{-1} μm^{-1}]'\n"; // TODO: add textual description (PDF?)
  gp << "set grid\n";
  gp << "set nokey\n";

  // FSSP range
  gp << "set arrow from .5," << ymin << " to .5," << ymax << " nohead\n";
  gp << "set arrow from 25," << ymin << " to 25," << ymax << " nohead\n";

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

      std::map<float, float> focus_d;
      std::map<float, float> focus_w;

      //info on the number and location of histogram edges
      vector<quantity<si::length>> left_edges_rd = bins_dry();
      int nsd = left_edges_rd.size() - 1;
      vector<quantity<si::length>> left_edges_rw = bins_wet();
      int nsw = left_edges_rw.size() - 1;

      for (int i = 0; i < nsd; ++i)
      {
	const string name = "rd_rng" + zeropad(i) + "_mom0";
	blitz::Array<float, 2> tmp_d(1e-6 * h5load(h5, name, at));

	focus_d[left_edges_rd[i] / 1e-6 / si::metres] = sum(tmp_d(
	  blitz::Range(x-1, x+1),
	  blitz::Range(y-1, y+1)
	)) 
	/ 9  // mean over 9 gridpoints
	/ ((left_edges_rd[i+1] - left_edges_rd[i]) / 1e-6 / si::metres); // per micrometre
      }

      for (int i = 0; i < nsw; ++i)
      {
	const string name = "rw_rng" + zeropad(i + off) + "_mom0";
	blitz::Array<float, 2> tmp_w(1e-6 * h5load(h5, name, at));

	focus_w[left_edges_rw[i] / 1e-6 / si::metres] = sum(tmp_w(
	  blitz::Range(x-1, x+1),
	  blitz::Range(y-1, y+1)
	)) 
	/ 9 
	/ ((left_edges_rw[i+1] - left_edges_rw[i]) / 1e-6 / si::metres); // per micrometre
      }

      notice_macro("setting-up plot parameters");
      gp << "plot"
	 << "'-' with histeps title 'wet radius' lw 3 lc rgb 'blue'," 
	 << "'-' with histeps title 'dry radius' lw 1 lc rgb 'red' " << endl;
      gp.send(focus_w);
      gp.send(focus_d);

      lbl -= 2;
    }
    lbl = 'j';
  }
}
