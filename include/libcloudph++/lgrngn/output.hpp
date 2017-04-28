#pragma once 
#include "chem.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    enum output_t // chem has to be first!
    {
      outHNO3 = HNO3, 
      outNH3 = NH3, 
      outCO2 = CO2, 
      outSO2 = SO2, 
      outH2O2 = H2O2, 
      outO3 = O3,
      outS_VI = S_VI,                        
      outH = H,                          
      outliq_vol,
      outdry_vol
    };
  };
};
