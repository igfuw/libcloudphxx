#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    enum chem_species_t 
    {
      HNO3, NH3, CO2, SO2, H2O2, O3,  // both gas and total dissolved chem species
                                      //    (must be first as gas species vectors use only these!)
      S_VI,                           // must be after above and before the below to use odeint 
                                      //    on a part of the vector only - odeint uses SO2 --> S_VI!
      H,                              // not used in the odeint vector  
                                             
      chem_gas_n   = O3 + 1,
      chem_rhs_beg = SO2,
      chem_rhs_fin = S_VI + 1,
      chem_all     = H + 1
    };
  };
};
