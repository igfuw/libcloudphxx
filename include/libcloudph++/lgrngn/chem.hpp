#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    enum chem_species_t
    {
      SO2, H2O2, O3, // both gas and aqueous (must be first as gas species vectors use only these!)
      H, OH, HSO3, SO3, S_VI, HSO4, SO4, // only aqueous
      chem_gas_n=O3+1,
      chem_aq_n=SO4+1
    };
  };
};
