#include "../common.hpp"
#include "../fig_a/gnuplot.hpp"
#include "../../src/icmw8_case1.hpp"

int main(int ac, char** av)
{
  std::map<std::string, int> n({{"x", 15},{"z", 15}});

  Gnuplot gp;
  init(gp, "./plot.svg", 1, 1, n);

  using namespace icmw8_case1;

  float dx = (X / si::metres) / (n["x"]-1), dz = (Z / si::metres) / (n["z"]-1);

  blitz::Array<float, 2> tmp(n["x"], n["z"]);

  float A = (w_max / si::metres_per_second) * (X / si::metres) / pi<float>();

  blitz::firstIndex ix;
  blitz::secondIndex jx;

  std::vector<std::vector<float>> v(4);

  tmp = (ix + .5) / n["x"] * (n["x"] - 1);

  for (const auto &i : tmp) v[0].push_back(i);
  
  tmp = (jx + .5) / n["z"] * (n["z"] - 1);

  for (const auto &i : tmp) v[1].push_back(i);

  tmp = - A * ( 
    psi((ix+.5)/n["x"], (jx+.5+.5)/n["z"])-
    psi((ix+.5)/n["x"], (jx+.5-.5)/n["z"])
  ) / dz                  // numerical derivative
  / rhod()((jx+.5) * dz); // psi defines rho_d times velocity

  for (const auto &i : tmp) v[2].push_back(i);

  tmp = A * ( 
    psi((ix+.5+.5)/n["x"], (jx+.5)/n["z"]) -
    psi((ix+.5-.5)/n["x"], (jx+.5)/n["z"])
  ) / dx  
  / rhod()(jx * dz);

  for (const auto &i : tmp) v[3].push_back(i);

  gp << "scl = .666\n";
  gp << "set title 'momentum field'\n";
  gp << "set xrange [0:" << n["x"]-1 << "]\n";
  gp << "set yrange [0:" << n["z"]-1 << "]\n";
  gp << "set cbrange [0:5]\n";
  gp << "splot '-' using 1:2:(0):(scl*$3):(scl*$4):(0):(sqrt($3**2 + $4**2)) with vectors linecolor rgbcolor variable notitle\n";
  gp.send(v);
}
