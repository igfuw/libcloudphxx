/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

// command-line option handling
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include "icmw8_case1.hpp"
namespace po = boost::program_options;

// some globals for option handling
int ac; 
char** av; // TODO: write it down to a file as in icicle ... write the default (i.e. not specified) values as well!
po::options_description opts_main("General options"); 

void handle_opts(
  po::options_description &opts_micro,
  po::variables_map &vm 
)
{
  opts_main.add(opts_micro);
  po::store(po::parse_command_line(ac, av, opts_main), vm); // could be exchanged with a config file parser

  // hendling the "help" option
  if (vm.count("help"))
  {
    std::cout << opts_main;
    exit(EXIT_SUCCESS);
  }
  po::notify(vm); // includes checks for required options
}

template<class real_t, class setup_t>
void setopts_common(
  setup_t &setup
)
{
  po::options_description opts("setup options");
  opts.add_options()
    ("th_0", po::value<real_t>()->default_value(289),  "th_0")
    ("rv_0", po::value<real_t>()->default_value(7.5e-3),  "th_0")
    ("p_0", po::value<real_t>()->default_value(101500),  "th_0")
    ("w_max", po::value<real_t>()->default_value(.6),  "th_0")
    ("z_0", po::value<real_t>()->default_value(0),  "th_0")
    ("Z", po::value<real_t>()->default_value(1500),  "th_0")
    ("X", po::value<real_t>()->default_value(1500),  "th_0")
    ("dt", po::value<real_t>()->default_value(1),  "th_0")
    ("mean_rd1", po::value<real_t>()->default_value(0.02e-6),  "th_0")
    ("mean_rd2", po::value<real_t>()->default_value(0.075e-6),  "th_0")
    ("sdev_rd1", po::value<real_t>()->default_value(1.4),  "th_0")
    ("sdev_rd2", po::value<real_t>()->default_value(1.6),  "th_0")
    ("n1_stp", po::value<real_t>()->default_value(60e6),  "th_0")
    ("n2_stp", po::value<real_t>()->default_value(40e6),  "th_0")
    ("kappa", po::value<real_t>()->default_value(.61),  "th_0")
    ("chem_b", po::value<real_t>()->default_value(.55),  "th_0")
    ("tau_rlx", po::value<real_t>()->default_value(300),  "th_0")
    ("z_rlx", po::value<real_t>()->default_value(200),  "th_0")
  ;
  po::variables_map vm;

  opts_main.add(opts);
  po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); // ignores unknown

  setup.th_0 = vm["th_0"].as<real_t>() * si::kelvins;
  setup.rv_0 = vm["rv_0"].as<real_t>();
  setup.p_0 = vm["p_0"].as<real_t>() * si::pascals;
  setup.w_max = vm["w_max"].as<real_t>() * si::metres_per_second;
  setup.z_0 = vm["z_0"].as<real_t>() * si::metres;
  setup.Z = vm["Z"].as<real_t>() * si::metres;
  setup.X = vm["X"].as<real_t>() * si::metres;
  setup.dt = vm["dt"].as<real_t>() * si::seconds;
  setup.mean_rd1 = vm["mean_rd1"].as<real_t>() * si::metres;
  setup.mean_rd2 = vm["mean_rd2"].as<real_t>() * si::metres;
  setup.sdev_rd1 = vm["sdev_rd1"].as<real_t>();
  setup.sdev_rd2 = vm["sdev_rd2"].as<real_t>();
  setup.n1_stp = vm["n1_stp"].as<real_t>() / si::cubic_metres;
  setup.n2_stp = vm["n2_stp"].as<real_t>() / si::cubic_metres;
  setup.kappa = vm["kappa"].as<real_t>();
  setup.chem_b = vm["chem_b"].as<real_t>();
  setup.tau_rlx = vm["tau_rlx"].as<real_t>() * si::seconds;
  setup.z_rlx = vm["z_rlx"].as<real_t>() * si::metres;
}





