#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
//<listing>
    enum class as_t { undefined, implicit, euler, pred_corr }; 
//</listing>
    const std::unordered_map<as_t, std::string> as_name = {
      {as_t::undefined, "undefined"},
      {as_t::implicit, "implicit"},
      {as_t::euler, "euler"},
      {as_t::pred_corr, "pred_corr"}
    };
  };
};
