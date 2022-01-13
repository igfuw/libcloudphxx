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
  };
};
