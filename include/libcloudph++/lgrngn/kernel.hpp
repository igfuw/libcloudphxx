#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace kernel_t //separate namespace to avoid member name conflicts with terminal_velocity enumerator, TODO: in c++11 change it to an enum class
    {
//<listing>
      enum kernel_t { undefined, geometric, golovin }; 
//</listing>
    };
  };
};
