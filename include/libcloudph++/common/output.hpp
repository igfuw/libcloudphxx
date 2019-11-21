#pragma once 
#include "chem.hpp"

namespace libcloudphxx
{
  namespace common
  {
    enum output_t // chem has to be first!
    {
      outHNO3 = chem::HNO3, 
      outNH3 = chem::NH3, 
      outCO2 = chem::CO2, 
      outSO2 = chem::SO2, 
      outH2O2 = chem::H2O2, 
      outO3 = chem::O3,
      outS_VI = chem::S_VI,                        
      outH = chem::H,                          
      outliq_vol,
      outdry_vol
    };

    const std::map<output_t, std::string> output_names
    {
      {outHNO3   , "HNO3"}, 
      {outNH3    , "NH3"}, 
      {outCO2    , "CO2"}, 
      {outSO2    , "SO2"}, 
      {outH2O2   , "H2O2"}, 
      {outO3     , "O3"},
      {outS_VI   , "S_VI"},                        
      {outH      , "H"},                          
      {outliq_vol, "liquid_volume"},
      {outdry_vol, "dry_volume"}
    };
  };
};
