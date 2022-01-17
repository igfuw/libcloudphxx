#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
//<listing>
    enum backend_t { undefined, serial, OpenMP, CUDA, multi_CUDA }; 
//</listing>
    const std::unordered_map<backend_t, std::string> backend_name = {
      {undefined, "undefined"},
      {serial, "serial"},
      {OpenMP, "OpenMP"},
      {CUDA, "CUDA"},
      {multi_CUDA, "multi_CUDA"}
    };
  };
};
