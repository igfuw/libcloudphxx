#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    enum chem_species_t
    {
      HNO3, NH3, CO2, SO2, H2O2, O3,         // both gas and dissolved chem species
                                             //    (must be first as gas species vectors use only these!)
      HSO3, S_VI, SO3,                       // must be after above and before the below to use odeint 
                                             //    on a part of the vector only - odeint uses SO2 --> SO3!
      H, OH, HSO4, SO4, HCO3, CO3, NH4, NO3, // not used in the odeint vector  
                                             
      chem_gas_n   = O3 + 1,
      chem_rhs_beg = SO2,
      chem_rhs_fin = SO3 + 1,
      chem_all     = NO3 + 1
    };
  };
};
