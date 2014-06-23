#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
// TODO: suffix with _t?
    enum chem_gas {
      gSO2, gO3, gH2O2,
      chem_gas_n=gH2O2+1
    };
    enum chem_aq 
    {
      H, OH, SO2, O3, H2O2, HSO3, SO3, S_VI, HSO4, SO4, 
      chem_aq_n=SO4+1
    };
  };
};
