#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/moist_air.hpp>

using namespace libcloudphxx::common::molar_mass;
using namespace libcloudphxx::common::moist_air;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string dir = string(av[1]) + "/thesis_2017/thesis_chap_6/";

  for (const std::string sim_run :
           {
            "out_case_base",
            "out_case1a",
            "out_case3",
            "out_case4",
            "out_case5",
            "out_case6"
           })
  {
    std::string h5  = dir + sim_run;
    auto n = h5n(h5);

    for (int at = 0; at < n["t"]; ++at)
    {
      if (at * n["outfreq"] >= 10000 || at == 0) //plot only after spinup
      {
        for (auto &plt : std::set<std::string>({"SO2g",  "O3g",    "H2O2g",  "CO2g",   "NH3g",  "HNO3g", 
                                               "S_IV_aq", "S_VI_aq", "O3_aq",  "H2O2_aq", "H_aq",  
                                               "C_IV_aq", "N_III_aq", "N_V_aq"}))

        {
          Gnuplot gp;
          init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".svg", 1, 1, n); 

          if (plt == "SO2g")
          {
            auto chem = h5load(h5, plt, at * n["outfreq"]) / (M_SO2<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'SO_{2} conc. [mol/μg of dry air]'\n";
            gp << "set cbrange [2.5 : 7.5]\n";
            plot(gp, chem);
          }
          else if (plt == "O3g")
          {
            auto chem = h5load(h5, plt, at * n["outfreq"]) / (M_O3<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'trace gas O_{3} conc. [mol/μg of dry air]'\n";
            gp << "set cbrange [871.8 : 873.4]\n";
            plot(gp, chem);
          }
          else if (plt == "H2O2g") 
          {
            auto chem = h5load(h5, plt, at * n["outfreq"]) / (M_H2O2<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'trace gas H_{2}O_{2} conc. [mol/μg of dry air]'\n";
            gp << "set cbrange [0 : 18]\n";
            plot(gp, chem);
          }
          else if (plt == "CO2g")
          {
            auto chem = h5load(h5, plt, at * n["outfreq"]) / (M_CO2<float>() * si::moles / si::kilograms) * 1e6;
            gp << "set title 'trace gas CO_{2} conc. [mol/mg of dry air]'\n";
            //gp << "set cbrange [360:361.2]\n";
            plot(gp, chem);
          }
          else if (plt == "NH3g") // 100e-12
          {
            auto chem = h5load(h5, plt, at * n["outfreq"]) / (M_NH3<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'trace gas NH_{3} conc. [mol/μg of dry air]'\n";
            //gp << "set cbrange [0:0.11]\n";
            plot(gp, chem);
          }
          else if (plt == "HNO3g") // 100e-12
          {
            auto chem = h5load(h5, plt, at * n["outfreq"]) / (M_HNO3<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'trace gas HNO_{3} conc. [mol/μg of dry air]'\n";
            //gp << "set cbrange [0:0.12]\n";
            plot(gp, chem);
          }

          else if (plt == "S_IV_aq")
          {
            auto chem = h5load(h5, "chem_S_IV_aq", at * n["outfreq"]) / (M_SO2_H2O<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'dissolved S^{IV} conc. [mol/μg of dry air]'\n";
            //gp << "set cbrange [0:1e-5]\n";
            plot(gp, chem);
          }
          else if (plt == "S_VI_aq")
          {
            auto chem = h5load(h5, "chem_S_VI_aq", at * n["outfreq"]) / (M_H2SO4<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'S^{VI} conc. [mol/μg of dry air]'\n";
            gp << "set cbrange [0 : 12]\n";
            plot(gp, chem);
          }
          else if (plt == "O3_aq")
          {
            auto chem = h5load(h5, "chem_O3_aq", at * n["outfreq"]) / (M_O3<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'dissolved O_{3} conc. [mol/μg of dry air]'\n";
            gp << "set cbrange [0 : 0.0006]\n";
            plot(gp, chem);
          }
          else if (plt == "H2O2_aq")
          {
            auto chem = h5load(h5, "chem_H2O2_aq", at * n["outfreq"]) / (M_H2O2<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'dissolved H_{2}O_{2} conc. [mol/μg of dry air]'\n";
            gp << "set cbrange [0 : 25]\n";
            plot(gp, chem);
          }
          else if (plt == "H_aq")
          {
            auto chem = h5load(h5, "chem_H_aq", at * n["outfreq"]) / (M_H<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'H^{+} [mol/μg of dry air]'\n";
            //gp << "set cbrange [0:0.025]\n";
            plot(gp, chem);
          }
          else if (plt == "C_IV_aq")
          {
            auto chem = h5load(h5, "chem_C_IV_aq", at * n["outfreq"]) / (M_CO2_H2O<float>() * si::moles / si::kilograms) * 1e6;
            gp << "set title 'dissolved C^{IV} conc. [mol/mg of dry air]'\n";
            //gp << "set cbrange [0:1.6]\n";
            plot(gp, chem);
          }
          else if (plt == "N_III_aq")
          {
            auto chem = h5load(h5, "chem_N_III_aq", at * n["outfreq"]) / (M_NH3_H2O<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'dissolved N^{-III} conc. [mol/μg of dry air]'\n";
            //gp << "set cbrange [0:0.45]\n";
            plot(gp, chem);
          }
          else if (plt == "N_V_aq")
          {
            auto chem = h5load(h5, "chem_N_V_aq", at * n["outfreq"]) / (M_HNO3<float>() * si::moles / si::kilograms) * 1e9;
            gp << "set title 'dissolved N^{V} conc. [mol/μg of dry air]'\n";
            //gp << "set cbrange [0:0.6]\n";
            plot(gp, chem);
          }

          else assert(false);
        } // var loop
      } // if for plotting only after sinup
    } // time loop
  } // sim loop
} // main

