#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
//<listing>
    enum class RH_formula_t { pv_cc,      // RH = pv / pvs with pvs from the Clausius-Clapeyron equation
                              rv_cc,      // RH = rv / rvs with rvs from the Clausius-Clapeyron equation
                              pv_tet,     // RH = pv / pvs with pvs from the Tetens equation
                              rv_tet      // RH = rv / rvs with rvs from the Tetens equation
                            }; 
//</listing>
    const std::unordered_map<RH_formula_t, std::string> RH_formula_name = {
      {RH_formula_t::pv_cc, "pv_cc"},
      {RH_formula_t::rv_cc, "rv_cc"},
      {RH_formula_t::pv_tet, "pv_tet"},
      {RH_formula_t::rv_tet, "rv_tet"}
    };
  };
};
