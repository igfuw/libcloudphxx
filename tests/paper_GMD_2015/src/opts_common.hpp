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
    ("th_0", po::value<float>()->default_value(289),  "th_0")
  ;
  po::variables_map vm;

  opts_main.add(opts);
  po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); // ignores unknown

  setup.th_0 = vm["th_0"].as<float>() * si::kelvins;
}




