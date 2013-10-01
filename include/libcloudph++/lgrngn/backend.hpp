#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    // to make inclusion of Thrust not neccesarry here
    enum backend_t { serial, OpenMP, CUDA }; 
  };
};
