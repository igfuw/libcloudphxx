#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace RH_formula_t //separate namespace to avoid member name conflicts with kernel enumerator, TODO: in c++11 change it to an enum class
    {   
//<listing>
      enum RH_formula_t { pv_cc,      // RH = pv / pvs with pvs from the Clausius-Clapeyron equation
                          rv_cc,      // RH = rv / rvs with rvs from the Clausius-Clapeyron equation
                          pv_tet,     // RH = pv / pvs with pvs from the Tetens equation
                          rv_tet      // RH = rv / rvs with rvs from the Tetens equation
                        }; 
//</listing>
    }; 
  };
};
