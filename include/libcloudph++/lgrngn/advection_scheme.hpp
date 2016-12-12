#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace as_t //separate namespace to avoid member name conflicts with kernel enumerator, TODO: in c++11 change it to an enum class
    {   
//<listing>
      enum as_t { undefined, implicit, euler, pred_corr }; 
//</listing>
    }; 
  };
};
