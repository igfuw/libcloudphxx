#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    // to make inclusion of Thrust not neccesarry here
//<listing>
    enum backend_t { serial, OpenMP, CUDA, multi_CUDA }; 
//</listing>
  };
};
