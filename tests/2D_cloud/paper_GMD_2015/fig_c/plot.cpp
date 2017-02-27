#include "../common.hpp"
#include "../fig_a/gnuplot.hpp"
#include "../fig_a/hdf5.hpp"
#include "../../src/icmw8_case1.hpp"

int main(int ac, char** av)
{
  //read w_max, X and Z
  std::string filename = string(av[1]) + "/paper_GMD_2015/fig_a/out_lgrngn/const.h5";
  H5::H5File h5f(filename, H5F_ACC_RDONLY);
  auto h5_group = h5f.openGroup("/");
  auto h5_data  = h5f.openDataSet("G");
  float X, Z, w_max;
  {
    auto attr_x = h5_group.openAttribute("X");
    attr_x.read(attr_x.getDataType(), &X);
    auto attr_z = h5_group.openAttribute("Z");
    attr_z.read(attr_z.getDataType(), &Z);
    auto attr_w = h5_group.openAttribute("w_max");
    attr_w.read(attr_w.getDataType(), &w_max);
  }
  //read rhod
  blitz::Array<float, 2> rhod(h5load(filename, "G", 0, false));

  std::map<std::string, int> n({{"x", 15},{"z", 15}});

  Gnuplot gp;
  init(gp, string(av[1]) + "/paper_GMD_2015/fig_c/plot.svg", 1, 1, n);

  float dx = X / (n["x"]-1), dz = Z / (n["z"]-1);

  blitz::Array<float, 2> tmp(n["x"], n["z"]);

  float A = w_max * X / pi<float>();

  blitz::firstIndex ix;
  blitz::secondIndex jx;

  std::vector<std::vector<float>> v(4);

  tmp = (ix + .5) / n["x"] * (n["x"] - 1);

  for (const auto &i : tmp) v[0].push_back(i);
  
  tmp = (jx + .5) / n["z"] * (n["z"] - 1);

  for (const auto &i : tmp) v[1].push_back(i);

  tmp = - A * ( 
    config::psi((ix+.5)/n["x"], (jx+.5+.5)/n["z"])-
    config::psi((ix+.5)/n["x"], (jx+.5-.5)/n["z"])
  ) / dz                  // numerical derivative
  / rhod(ix, jx);         // psi defines rho_d times velocity
//  / rhod()((jx+.5) * dz);  //TODO - should interpolate to jx + .5

  for (const auto &i : tmp) v[2].push_back(i);

  tmp = A * ( 
    config::psi((ix+.5+.5)/n["x"], (jx+.5)/n["z"]) -
    config::psi((ix+.5-.5)/n["x"], (jx+.5)/n["z"])
  ) / dx  
  / rhod(ix, jx);
//  / rhod()(jx * dz);

  for (const auto &i : tmp) v[3].push_back(i);

  gp << "scl = .666\n";
  gp << "set title 'momentum field'\n";
  gp << "set xrange [0:" << n["x"]-1 << "]\n";
  gp << "set yrange [0:" << n["z"]-1 << "]\n";
  gp << "set cbrange [0:5]\n";
  gp << "splot '-' using 1:2:(0):(scl*$3):(scl*$4):(0):(sqrt($3**2 + $4**2)) with vectors linecolor rgbcolor variable notitle\n";
  gp.send(v);

}
