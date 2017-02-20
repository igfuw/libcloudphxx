/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include <boost/assign/ptr_map_inserter.hpp>  // for 'ptr_map_insert()'

#include "opts_common.hpp"
#include "kin_cloud_2d_lgrngn.hpp"
#include "kin_cloud_2d_lgrngn_chem.hpp"

// string parsing
#include <boost/spirit/include/qi.hpp>    
#include <boost/fusion/adapted/std_pair.hpp> 
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

/*
// a riddle for tired mind :) - how to make it work without passing ct_params

template <class solver_tp>
struct setopts_micro_chem_f;

template <template<class> class solver_tp, class ct_params_t>
struct setopts_micro_chem_f<solver_tp<ct_params_t>>
{
   void operator()(typename solver_tp<ct_params_t>::rt_params_t &rt_params,
                   typename std::enable_if<std::is_same<solver_tp<ct_params_t>, kin_cloud_2d_lgrngn<ct_params_t>>::value>::type* = 0                  )
   {
     using solver_t = solver_tp<ct_params_t>;
     rt_params.outvars = {
       {solver_t::ix::th, {"th", "[K]"}},
       {solver_t::ix::rv, {"rv", "[kg kg-1]"}}
     };
   }
   
   void operator()(typename solver_tp<ct_params_t>::rt_params_t &rt_params,
                   typename std::enable_if<std::is_same<solver_tp<ct_params_t>, kin_cloud_2d_lgrngn_chem<ct_params_t>>::value>::type* = 0                  )
   {
     using solver_t = solver_tp<ct_params_t>;
     rt_params.outvars = {
       {solver_t::ix::th, {"th", "[K]"}},
       {solver_t::ix::rv, {"rv", "[kg kg-1]"}},
       {solver_t::ix::SO2g,  {"SO2g", "[dimensionless]"}},
       {solver_t::ix::O3g,   {"O3g",  "[dimensionless]"}},
       {solver_t::ix::H2O2g, {"H2O2g","[dimensionless]"}},
       {solver_t::ix::CO2g,  {"CO2g", "[dimensionless]"}},
       {solver_t::ix::NH3g,  {"NH3g", "[dimesnionless]"}},
       {solver_t::ix::HNO3g, {"HNO3g","[dimensionless]"}}
     };
   }
};

//and then to call it ...
  setopts_micro_chem_f<solver_t> chem_f{};
  chem_f(rt_params);

//increment the number of hours you have spent on it: 1
*/
template <class solver_t, class ct_params_t>
void parse_moms(
  std::string opt, 
  po::variables_map vm, 
  typename solver_t::rt_params_t &rt_params,
  typename std::enable_if<std::is_same<solver_t, kin_cloud_2d_lgrngn<ct_params_t>>::value>::type* = 0
) {
  namespace qi = boost::spirit::qi;
  namespace phoenix = boost::phoenix;

  std::string val = vm[opt].as<std::string>();
  auto first = val.begin();
  auto last  = val.end();

  std::vector<std::pair<std::string, std::string>> min_maxnum;
  outmom_t<config::real_t> &moms =
    opt == "out_dry" ? rt_params.out_dry :
      rt_params.out_wet;

  const bool result = qi::phrase_parse(first, last, 
    *(
      *(qi::char_-":")  >>  qi::lit(":") >>  
      *(qi::char_-";")  >> -qi::lit(";") 
    ),
    boost::spirit::ascii::space, min_maxnum
  );    
  if (!result || first != last) BOOST_THROW_EXCEPTION(po::validation_error(
      po::validation_error::invalid_option_value, opt, val 
  ));  

  for (auto &ss : min_maxnum)
  {
    int sep = ss.second.find('|'); 

    moms.push_back(outmom_t<config::real_t>::value_type({
      outmom_t<config::real_t>::value_type::first_type(
        boost::lexical_cast<config::real_t>(ss.first) * si::metres,
        boost::lexical_cast<config::real_t>(ss.second.substr(0, sep)) * si::metres
      ), 
      outmom_t<config::real_t>::value_type::second_type()
    }));

    // TODO catch (boost::bad_lexical_cast &)

    std::string nums = ss.second.substr(sep+1);;
    auto nums_first = nums.begin();
    auto nums_last  = nums.end();

    const bool result = qi::phrase_parse(
      nums_first, 
      nums_last, 
      (
        qi::int_[phoenix::push_back(phoenix::ref(moms.back().second), qi::_1)]
            >> *(',' >> qi::int_[phoenix::push_back(phoenix::ref(moms.back().second), qi::_1)])
      ),
      boost::spirit::ascii::space
    );    
    if (!result || nums_first != nums_last) BOOST_THROW_EXCEPTION(po::validation_error(
        po::validation_error::invalid_option_value, opt, val // TODO: report only the relevant part?
    ));  
  }
}
template <class solver_t, class ct_params_t>
void parse_moms(
  std::string opt, 
  po::variables_map vm, 
  typename solver_t::rt_params_t &rt_params,
  typename std::enable_if<std::is_same<solver_t, kin_cloud_2d_lgrngn_chem<ct_params_t>>::value>::type* = 0
) {
  namespace qi = boost::spirit::qi;
  namespace phoenix = boost::phoenix;

  std::string val = vm[opt].as<std::string>();

  auto first = val.begin();
  auto last  = val.end();

  std::vector<std::pair<std::string, std::string>> min_maxnum;

  outmom_t<config::real_t> &moms =
    opt == "out_dry" ? rt_params.out_dry :
      opt == "out_wet" ? rt_params.out_wet : 
         opt == "out_chem" ? rt_params.out_chem :
           rt_params.out_wet_pH;

  const bool result = qi::phrase_parse(first, last, 
    *(
      *(qi::char_-":")  >>  qi::lit(":") >>  
      *(qi::char_-";")  >> -qi::lit(";") 
    ),
    boost::spirit::ascii::space, min_maxnum
  );    
  if (!result || first != last) BOOST_THROW_EXCEPTION(po::validation_error(
      po::validation_error::invalid_option_value, opt, val 
  ));  

  for (auto &ss : min_maxnum)
  {
    int sep = ss.second.find('|'); 

    moms.push_back(outmom_t<config::real_t>::value_type({
      outmom_t<config::real_t>::value_type::first_type(
        boost::lexical_cast<config::real_t>(ss.first) * si::metres,
        boost::lexical_cast<config::real_t>(ss.second.substr(0, sep)) * si::metres
      ), 
      outmom_t<config::real_t>::value_type::second_type()
    }));

    // TODO catch (boost::bad_lexical_cast &)

    std::string nums = ss.second.substr(sep+1);;
    auto nums_first = nums.begin();
    auto nums_last  = nums.end();

    const bool result = qi::phrase_parse(
      nums_first, 
      nums_last, 
      (
        qi::int_[phoenix::push_back(phoenix::ref(moms.back().second), qi::_1)]
            >> *(',' >> qi::int_[phoenix::push_back(phoenix::ref(moms.back().second), qi::_1)])
      ),
      boost::spirit::ascii::space
    );    
    if (!result || nums_first != nums_last) BOOST_THROW_EXCEPTION(po::validation_error(
        po::validation_error::invalid_option_value, opt, val // TODO: report only the relevant part?
    ));  
  }
}

template <class solver_t, class ct_params_t>
void setopts_micro_chem(
  typename solver_t::rt_params_t &rt_params, 
  std::set<std::string> &out_set,
  typename std::enable_if<std::is_same<solver_t, kin_cloud_2d_lgrngn<ct_params_t>>::value>::type* = 0
) {//TODO - add out_wet and others here
  rt_params.outvars = {
    {solver_t::ix::th, {"th", "[K]"}},
    {solver_t::ix::rv, {"rv", "[kg kg-1]"}}
  };
  out_set = {"out_dry", "out_wet"};
}

template <class solver_t, class ct_params_t>
void setopts_micro_chem(
  typename solver_t::rt_params_t &rt_params, 
  std::set<std::string> &out_set,
  typename std::enable_if<std::is_same<solver_t, kin_cloud_2d_lgrngn_chem<ct_params_t>>::value>::type* = 0
) {
  rt_params.outvars = {
    {solver_t::ix::th, {"th", "[K]"}},
    {solver_t::ix::rv, {"rv", "[kg kg-1]"}},
    {solver_t::ix::SO2g,  {"SO2g", "[dimensionless]"}},
    {solver_t::ix::O3g,   {"O3g",  "[dimensionless]"}},
    {solver_t::ix::H2O2g, {"H2O2g","[dimensionless]"}},
    {solver_t::ix::CO2g,  {"CO2g", "[dimensionless]"}},
    {solver_t::ix::NH3g,  {"NH3g", "[dimesnionless]"}},
    {solver_t::ix::HNO3g, {"HNO3g","[dimensionless]"}}
  };
  out_set = {"out_dry", "out_wet", "out_chem", "out_wet_pH"};
}

// simulation and output parameters for micro=lgrngn
template <class solver_t, class ct_params_t>
void setopts_micro(
  typename solver_t::rt_params_t &rt_params, 
  int nx, int nz, int nt, config::setup_t &setup,
  typename std::enable_if<std::is_same<
    decltype(solver_t::rt_params_t::cloudph_opts),
    libcloudphxx::lgrngn::opts_t<typename solver_t::real_t>
  >::value>::type* = 0
)
{
  po::options_description opts("Lagrangian microphysics options"); 
  opts.add_options()
    ("backend", po::value<std::string>()->required() , "one of: CUDA, OpenMP, serial")
    ("async", po::value<bool>()->default_value(true), "use CPU for advection while GPU does micro (ignored if backend != CUDA)")
    ("sd_conc", po::value<unsigned long long>()->required() , "super-droplet number per grid cell (unsigned long long)")
    // processes
    ("adve", po::value<bool>()->default_value(rt_params.cloudph_opts.adve) , "particle advection     (1=on, 0=off)")
    ("sedi", po::value<bool>()->default_value(rt_params.cloudph_opts.sedi) , "particle sedimentation (1=on, 0=off)")
    ("cond", po::value<bool>()->default_value(rt_params.cloudph_opts.cond) , "condensational growth  (1=on, 0=off)")
    ("coal", po::value<bool>()->default_value(rt_params.cloudph_opts.coal) , "collisional growth     (1=on, 0=off)")
    ("rcyc", po::value<bool>()->default_value(rt_params.cloudph_opts.rcyc) , "recycling of droplets  (1=on, 0=off)")
    ("chem_dsl", po::value<bool>()->default_value(rt_params.cloudph_opts.chem_dsl) , "dissolving trace gases (1=on, 0=off)")
    ("chem_dsc", po::value<bool>()->default_value(rt_params.cloudph_opts.chem_dsc) , "dissociation           (1=on, 0=off)")
    ("chem_rct", po::value<bool>()->default_value(rt_params.cloudph_opts.chem_rct) , "aqueous chemistry      (1=on, 0=off)")
    ("chem_switch", po::value<bool>()->default_value(rt_params.cloudph_opts_init.chem_switch) , "aqueous chemistry (1=on, 0=off)")
    // free parameters
    ("sstp_cond", po::value<int>()->default_value(rt_params.cloudph_opts_init.sstp_cond), "no. of substeps for condensation")
    ("sstp_coal", po::value<int>()->default_value(rt_params.cloudph_opts_init.sstp_coal), "no. of substeps for coalescence")
    ("sstp_chem", po::value<int>()->default_value(rt_params.cloudph_opts_init.sstp_chem), "no. of substeps for chemistry")
    ("dev_count", po::value<int>()->default_value(rt_params.cloudph_opts_init.dev_count), "no of GPUs to use")
    ("rng_seed",  po::value<int>()->default_value(rt_params.cloudph_opts_init.rng_seed), "seed for random super droplet init")
    // output
    ("out_dry", po::value<std::string>()->default_value("0:1|0"),       "dry radius ranges and moment numbers (r1:r2|n1,n2...;...)")
    ("out_wet", po::value<std::string>()->default_value(".5e-6:25e-6|0,1,2,3;25e-6:1|0,3,6"),  "wet radius ranges and moment numbers (r1:r2|n1,n2...;...)")
    ("out_wet_pH", po::value<std::string>()->default_value("0:1|0"), "wet radius ranges for output of H+ and S_VI)")
    ("out_chem", po::value<std::string>()->default_value("0:1|0"),   "dry radius ranges for which chem mass is outputted")
    // collision and sedimentation
    ("kernel", po::value<std::string>()->default_value("geometric"), "collision kernel (geometric, long, hall, hall_davis_no_waals, golovin, onishi_hall, onishi_hall_davis_no_waals, vohl_davis_no_waals, hall_pinsky_cumulonimbus, hall_pinsky_stratocumulus)")
    ("terminal_velocity", po::value<std::string>()->default_value("khvorostyanov_spherical"), "sedimentation velocity (khvorostyanov_spherical, khvorostyanov_nonspherical, beard76, beard77, beard77fast)")
    ("adve_scheme", po::value<std::string>()->default_value("implicit"), "advection for super-droplets (implicit, euler, pred_corr)")
    // TODO: MAC, HAC, vent_coef
  ;
  po::variables_map vm;
  handle_opts(opts, vm);
      
  std::string backend_str = vm["backend"].as<std::string>();
  if (backend_str == "CUDA") rt_params.backend = libcloudphxx::lgrngn::CUDA;
  else if (backend_str == "OpenMP") rt_params.backend = libcloudphxx::lgrngn::OpenMP;
  else if (backend_str == "serial") rt_params.backend = libcloudphxx::lgrngn::serial;
  else if (backend_str == "multi_CUDA") rt_params.backend = libcloudphxx::lgrngn::multi_CUDA;

  rt_params.async = vm["async"].as<bool>();

  rt_params.cloudph_opts_init.sd_conc = vm["sd_conc"].as<unsigned long long>();
  rt_params.cloudph_opts_init.nx = nx;
  rt_params.cloudph_opts_init.nz = nz;
  if (backend_str == "multi_CUDA")
    rt_params.cloudph_opts_init.n_sd_max = 1.5 * nx *  nz * rt_params.cloudph_opts_init.sd_conc;
  else
    rt_params.cloudph_opts_init.n_sd_max = nx *  nz * rt_params.cloudph_opts_init.sd_conc;
 
  boost::assign::ptr_map_insert<
    config::log_dry_radii<config::real_t> // value type
  >(
    rt_params.cloudph_opts_init.dry_distros // map
  )(
    setup.kappa, // key
    setup
  );

  // process toggling
  rt_params.cloudph_opts.adve = vm["adve"].as<bool>();
  rt_params.cloudph_opts.sedi = vm["sedi"].as<bool>();
  rt_params.cloudph_opts.cond = vm["cond"].as<bool>();
  rt_params.cloudph_opts.coal = vm["coal"].as<bool>();

  rt_params.cloudph_opts.rcyc = vm["rcyc"].as<bool>();
  rt_params.cloudph_opts.chem_dsl = vm["chem_dsl"].as<bool>();
  rt_params.cloudph_opts.chem_dsc = vm["chem_dsc"].as<bool>();
  rt_params.cloudph_opts.chem_rct = vm["chem_rct"].as<bool>();
  rt_params.cloudph_opts_init.chem_switch = vm["chem_switch"].as<bool>();
  rt_params.cloudph_opts_init.chem_rho    = setup.chem_rho * si::cubic_metres / si::kilograms;

  // free parameters
  rt_params.cloudph_opts_init.sstp_cond = vm["sstp_cond"].as<int>();
  rt_params.cloudph_opts_init.sstp_coal = vm["sstp_coal"].as<int>();
  rt_params.cloudph_opts_init.sstp_chem = vm["sstp_chem"].as<int>();
  rt_params.cloudph_opts_init.dev_count = vm["dev_count"].as<int>();
  rt_params.cloudph_opts_init.rng_seed  = vm["rng_seed"].as<int>();

  // advection of super droplets choice
  if (vm["adve_scheme"].as<std::string>() == "implicit") {
    rt_params.cloudph_opts_init.adve_scheme = libcloudphxx::lgrngn::as_t::implicit;
  }
  if (vm["adve_scheme"].as<std::string>() == "euler") {
    rt_params.cloudph_opts_init.adve_scheme = libcloudphxx::lgrngn::as_t::euler;
  }
  if (vm["adve_scheme"].as<std::string>() == "pred_corr") {
    rt_params.cloudph_opts_init.adve_scheme = libcloudphxx::lgrngn::as_t::pred_corr;
  }
   else {
    std::cerr<<"Invalid advection choice"<<std::endl;
    assert(false); 
  }

  // coalescence kernel choice
  if (vm["kernel"].as<std::string>() == "geometric") {
    rt_params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::geometric;
    // halving the collection efficiency to match the timing of precipitation onset in the blk_2m scheme
    rt_params.cloudph_opts_init.kernel_parameters = {.5}; 
  }
  else if (vm["kernel"].as<std::string>() == "long") {
    // akin to geometric but a bit more sophisticated  (not realistic)
    rt_params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::Long;
  }
  else if (vm["kernel"].as<std::string>() == "hall") {
    // geometric kernel + collision efficiencies taken from lookup table (bad lookups for droplets smaller than 10 microns)
    rt_params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::hall;
  }
  else if (vm["kernel"].as<std::string>() == "hall_davis_no_waals") { 
    // geometric kernel + collision efficiencies taken from lookup table (like hall but collision efficiencies for small droplets are taken from Davis)
    rt_params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::hall_davis_no_waals;
  }
  else if (vm["kernel"].as<std::string>() == "golovin") {
    // a kernel that allows for analictic solution for collisons (good for testing, not for cloud simulations)
    rt_params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::golovin;
    rt_params.cloudph_opts_init.kernel_parameters = {1.}; // efficiency parameter 
  }
  else if (vm["kernel"].as<std::string>() == "onishi_hall") {
    // like hall but also taking into acount turbulent effects (modifying kernel and collision efficiencies) 
    // it's not geometric kernel anymore
    rt_params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::onishi_hall;
    rt_params.cloudph_opts_init.kernel_parameters = {0.01, 100}; // (rate of dissipation epsilon [m^2/s^3], Taylor microscale Reynolds number)
                                                  // 0.04 - cumulus
                                                  // 0.01 - stratocumulus
  }
  else if (vm["kernel"].as<std::string>() == "onishi_hall_davis_no_waals") {
    // like onisi_hall but collision efficiencies taken from Davis
    rt_params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::onishi_hall_davis_no_waals;
    rt_params.cloudph_opts_init.kernel_parameters = {0.01, 100}; // (rate of dissipation epsilon [m^2/s^3], Taylor microscale Reynolds number)
  }
  else if (vm["kernel"].as<std::string>() == "vohl_davis_no_waals") {
    // geometric kernel collision efficiencies for big droplets are taken from Vohl and for small from Davis
    rt_params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::vohl_davis_no_waals;
  }
  else if (vm["kernel"].as<std::string>() == "hall_pinsky_cumulonimbus") {
    // like geometric, collision efficiencies for big droplets from Hall for small droplets from Pinski (taking into account turbulence in cumulonimbus) 
    rt_params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::hall_pinsky_cumulonimbus;
  }
  else if (vm["kernel"].as<std::string>() == "hall_pinsky_stratocumulus") {
    // like geometric, collision efficiencies for big droplets from Hall for small droplets from Pinski (taking into account turbulence in cumulonimbus) 
    rt_params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::hall_pinsky_stratocumulus;
  }
  else {
    std::cerr<<"Invalid kernel choice"<<std::endl;
    assert(false); 
  }

  // terminal velocity choice
  if (vm["terminal_velocity"].as<std::string>() == "khvorostyanov_spherical") {
    // took results for Beard and provided different fit (this is good only for small droplets)
    rt_params.cloudph_opts_init.terminal_velocity = libcloudphxx::lgrngn::vt_t::khvorostyanov_spherical;
  }
  else if (vm["terminal_velocity"].as<std::string>() == "khvorostyanov_nonspherical") {
    // took results for Beard and provided different fit (good only for big droplets)
    rt_params.cloudph_opts_init.terminal_velocity = libcloudphxx::lgrngn::vt_t::khvorostyanov_nonspherical;
  }
  else if (vm["terminal_velocity"].as<std::string>() == "beard76") {
    // fit to measurements (complicated equation)
    rt_params.cloudph_opts_init.terminal_velocity = libcloudphxx::lgrngn::vt_t::beard76;
  } 
  else if (vm["terminal_velocity"].as<std::string>() == "beard77") {
    // simplified fit to Beard76
    rt_params.cloudph_opts_init.terminal_velocity = libcloudphxx::lgrngn::vt_t::beard77;
  } 
  else if (vm["terminal_velocity"].as<std::string>() == "beard77fast") {
    // creates a lookup table for p=1000hPa and scales the results for other values of pressure with some simple equation
    rt_params.cloudph_opts_init.terminal_velocity = libcloudphxx::lgrngn::vt_t::beard77fast;
  } 
  else {
    std::cerr<<"Invalid terminal velocity choice"<<std::endl;
    assert(false);
  }

  std::set<std::string> out_set;
  setopts_micro_chem<solver_t, ct_params_t>(rt_params, out_set);

  // parsing --out_dry and --out_wet options values
  // the format is: "rmin:rmax|0,1,2;rmin:rmax|3;..."
  for (auto &opt : out_set)
  {
    parse_moms<solver_t, ct_params_t>(opt, vm, rt_params);
  }
}
