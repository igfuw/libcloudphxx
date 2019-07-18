#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace SGS_length_scale_t //separate namespace to avoid member name conflicts with kernel enumerator, TODO: in c++11 change it to an enum class
    {   
//<listing>
      enum SGS_length_scale_t { vertical, geometric_mean, arithmetic_mean }; 
//</listing>
    }; 
  };
};
