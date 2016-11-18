#include <cstdlib> // system()
#include <set>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"

using std::ostringstream;
using std::set;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string opts_common = 
    "--outfreq=200 --nt=9000 --spinup=7200 --nx=76 --nz=76 --relax_th_rv=false";
  set<string> opts_micro({
    "--micro=blk_1m --outdir=out_blk_1m",
    "--micro=blk_2m --outdir=out_blk_2m"
  });

  for (auto &opts_m : opts_micro)
  {
    ostringstream cmd;
    cmd << av[1] << "/src/icicle " << opts_common << " " << opts_m;  
    notice_macro("about to call: " << cmd.str())

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
  }
}
