#pragma once

#include <blitz/array.h>
#include <gnuplot-iostream.h>
#include <map>

void init(
  Gnuplot &gp, 
  const std::string &file, 
  const int &ny, const int &nx, 
  std::map<std::string, int> n
)
{
  boost::filesystem::create_directories(
    boost::filesystem::path(file).parent_path()
  );

  gp << "set term svg dynamic enhanced fsize 14 size " << nx * 500 << "," << ny * 500 << "\n";
  gp << "set size square\n";
  gp << "set encoding utf8\n";
  // progressive-rock connoisseur palette ;)
  // gp << "set palette defined (0 '#FFFFFF', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')\n";
  //
  // Python matplotlib viridis palette 
  //   see https://github.com/BIDS/colormap/blob/master/colormaps.py for discussion
  //   and thank Nathaniel J. Smith, Stefan van der Walt and Eric Firing for the best palette ever.
  // Values for gnuplot were taken from https://github.com/Gnuplotting/gnuplot-palettes
  //   thanks are due to Hagen Wierstorf.
  // Palette was modified to start with white and not yellow
  gp << "set palette defined (0 '#ffffff', 1 '#aadc32', 2 '#5cc863', 3 '#27ad81', 4 '#21908d', 5 '#2c718e', 6 '#3b518b', 7 '#472c7a', 8 '#440154')\n";
  gp << "set view map\n";
  gp << "dx = 1500./" << n["x"]-1 << "\n"; 
  gp << "dy = 1500./" << n["z"]-1 << "\n"; 
  gp << "set xtics out scale .5 rotate by 60 ('0.0' 0, '0.3' 300/dx, '0.6' 600/dx, '0.9' 900/dx, '1.2' 1200/dx, '1.5' 1500/dx)\n"; 
  gp << "set xlabel 'x [km]'\n";
  gp << "set ytics out scale .5 rotate by 60 ('0.0' 0, '0.3' 300/dy, '0.6' 600/dy, '0.9' 900/dy, '1.2' 1200/dy, '1.5' 1500/dy)\n"; 
  gp << "set ylabel 'z [km]'\n";
  gp << "set output '" << file << "'\n";
  gp << "set grid\n";
  gp << "set multiplot layout " << ny << "," << nx << "\n";
}

template <class data_t>
void plot(Gnuplot &gp, const data_t &data)
{
  blitz::Array<float, 2> tmp(data);

  gp << "set xrange [0:" << tmp.extent(0)-1 << "]\n";
  gp << "set yrange [0:" << tmp.extent(1)-1 << "]\n";
  gp << "splot '-' binary" << gp.binfmt(tmp) << " scan=yx origin=(0,0,0) with image failsafe notitle\n";

  gp.sendBinary(tmp);
}
