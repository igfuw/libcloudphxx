#include <cstdlib> // system()
#include <list>
#include <map>
#include <string>
#include <sstream> // std::ostringstream
#include <fstream> 
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>

#include "../common.hpp"

using std::ostringstream;
using std::list;
using std::string;
using std::map;
using std::pair;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  int nt_load = 1, nt_calc = 5, n_rpt = 3;

  using str = std::string;
  using mss = std::map<str,str>;

  map<str,mss> proc({
    pair<str,mss>({"blk_1m", mss({
      pair<str,str>({"a___","--cond=off --cevp=off --revp=off --conv=off --accr=off --sedi=off"}),
      pair<str,str>({"ac__","--cond=on  --cevp=on  --revp=on  --conv=off --accr=off --sedi=off"}),
      pair<str,str>({"acc_","--cond=on  --cevp=on  --revp=on  --conv=on  --accr=on  --sedi=off"}),
      pair<str,str>({"accs","--cond=on  --cevp=on  --revp=on  --conv=on  --accr=on  --sedi=on "})
    })}),
    pair<str,mss>({"blk_2m", mss({
      pair<str,str>({"a___","--acti=off --cond=off --accr=off --acnv=off --sedi=off"}),
      pair<str,str>({"ac__","--acti=on  --cond=on  --accr=off --acnv=off --sedi=off"}),
      pair<str,str>({"acc_","--acti=on  --cond=on  --accr=on  --acnv=on  --sedi=off"}),
      pair<str,str>({"accs","--acti=on  --cond=on  --accr=on  --acnv=on  --sedi=on "})
    })}),
    pair<str,mss>({"lgrngn", mss({
      pair<str,str>({"a___","--adve=on --cond=off --coal=off --sedi=off"}),
      pair<str,str>({"ac__","--adve=on --cond=on  --coal=off --sedi=off"}),
      pair<str,str>({"acc_","--adve=on --cond=on  --coal=on  --sedi=off"}),
      pair<str,str>({"accs","--adve=on --cond=on  --coal=on  --sedi=on "}),
    })})
  });

  // output file
  std::ofstream output("timings.txt");
  if (!output.good()) throw std::runtime_error("failed to open output file");

  for (auto &micro : list<string>({"blk_1m", "blk_2m", "lgrngn"}))
  {
    for (auto &backend : micro == "lgrngn"
      ? list<string>({
        "--backend=CUDA", 
        "--backend=OpenMP",
        "--backend=serial"
      }) 
      : list<string>({""})
    ) {
      for (auto &env : micro == "lgrngn" && backend == "--backend=OpenMP"
	? list<string>({
          "OMP_NUM_THREADS=2",
          "OMP_NUM_THREADS=4"
        }) 
	: list<string>({""})
      ) {
	for (auto &sd_conc : micro == "lgrngn"
	  ? list<string>({
	    "--sd_conc=8",
	    "--sd_conc=32",
	    "--sd_conc=128"
	  }) 
	  : list<string>({""})
	) {
          // gnuplot-understendable output
          output << "# " << micro << " " << backend << " " << env << " " << sd_conc << std::endl;

	  for (auto prcs : list<string>({"a___","ac__","acc_","accs"}))
	  {
	    // multiplynig the time of the simulation to avoid measuring too short wall times
	    int mlt = 1;
            if (micro != "lgrngn") mlt *= 10;
            if (sd_conc == "--sd_conc=8") mlt *= 6; 
            if (sd_conc == "--sd_conc=32") mlt *= 3;
            if (backend == "--backend=CUDA") mlt *= 3;

	    double time_avg = 0;
	    for (auto &nt : list<int>({nt_load, nt_calc}))
	    {
	      ostringstream cmd;
	      cmd 
		<< env << " "
		<< av[1] 
		<< "/src/icicle --adv_serial=true --outdir=" 
                  << boost::filesystem::temp_directory_path().string() 
                  << "/" << boost::filesystem::unique_path("icicle-fig_b-%%%%%%%%").string() 
		<< " --outfreq=" << 44 * nt * mlt << " --nt=" << nt * mlt
		<< " --micro=" << micro << " " << backend << " " << sd_conc << " " << proc.at(micro).at(prcs);  

	      if (micro == "lgrngn") cmd << " --sstp_cond=10 --sstp_coal=10 --spinup=0"; // 44

	      notice_macro("about to call: " << cmd.str())

	      double time_min = std::numeric_limits<double>::max();
	      boost::timer::cpu_timer tmr;
	      for (int rpt=0; rpt < n_rpt; ++rpt) 
	      {
		tmr.start();
		if (EXIT_SUCCESS != system(cmd.str().c_str())) 
		  error_macro("model run failed: " << cmd.str())
		tmr.stop();
		double time = double(tmr.elapsed().wall) * 1e-9;
		if (time < time_min) time_min = time;
	      }
	      time_avg += (nt == nt_load ? -1 : 1) * time_min;
	    }
	    time_avg /= mlt * (nt_calc - nt_load);

            // gnuplot-understendable output
            output << time_avg << std::endl;
            output.flush();
	  }

          // gnuplot-understendable output
          output << std::endl;
          output << std::endl;
	}
      }
    }
  }
}
