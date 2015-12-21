/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "opts_common.hpp"
#include "kin_cloud_2d_blk_1m.hpp"

// simulation and output parameters for micro=blk_1m
template <class solver_t>
void setopts_micro(
  typename solver_t::rt_params_t &rt_params, 
  int nx, int nz, int nt,
  typename std::enable_if<std::is_same<
    decltype(solver_t::rt_params_t::cloudph_opts),
    libcloudphxx::blk_1m::opts_t<typename solver_t::real_t>
  >::value>::type* = 0
)
{
  po::options_description opts("Single-moment bulk microphysics options"); 
  opts.add_options()
    ("cond", po::value<bool>()->default_value(rt_params.cloudph_opts.cond) , "cloud water condensation (1=on, 0=off)")
    ("cevp", po::value<bool>()->default_value(rt_params.cloudph_opts.cevp) , "cloud water evaporation (1=on, 0=off)")
    ("revp", po::value<bool>()->default_value(rt_params.cloudph_opts.revp) , "rain water evaporation (1=on, 0=off)")
    ("conv", po::value<bool>()->default_value(rt_params.cloudph_opts.conv) , "autoconversion of cloud water into rain (1=on, 0=off)")
    ("accr", po::value<bool>()->default_value(rt_params.cloudph_opts.accr) , "cloud water collection by rain (1=on, 0=off)")
    ("sedi", po::value<bool>()->default_value(rt_params.cloudph_opts.sedi) , "rain water sedimentation (1=on, 0=off)")
//TODO: autoconv_threshold, epsilon
  ;
  po::variables_map vm;
  handle_opts(opts, vm);

  // Kessler scheme options
  rt_params.cloudph_opts.cond = vm["cond"].as<bool>();
  rt_params.cloudph_opts.cevp = vm["cevp"].as<bool>();
  rt_params.cloudph_opts.revp = vm["revp"].as<bool>();
  rt_params.cloudph_opts.conv = vm["conv"].as<bool>();
  rt_params.cloudph_opts.accr = vm["accr"].as<bool>();
  rt_params.cloudph_opts.sedi = vm["sedi"].as<bool>();

  // output variables
  rt_params.outvars = {
    // <TODO>: make it common among all three micro?
    {solver_t::ix::th, {"th", "[K]"}},
    {solver_t::ix::rv, {"rv", "[kg kg-1]"}},
    // </TODO>
    {solver_t::ix::rc, {"rc", "[kg kg-1]"}},
    {solver_t::ix::rr, {"rr", "[kg kg-1]"}}
  };
}
