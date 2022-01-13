#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
//<listing>
    enum class vt_t { undefined, beard76, beard77, beard77fast, khvorostyanov_spherical, khvorostyanov_nonspherical }; 
//</listing>
    const std::unordered_map<vt_t, std::string> vt_name = {
      {vt_t::undefined, "undefined"},
      {vt_t::beard76, "beard76"},
      {vt_t::beard77, "beard77"},
      {vt_t::beard77fast, "beard77fast"},
      {vt_t::khvorostyanov_spherical, "khvorostyanov_spherical"},
      {vt_t::khvorostyanov_nonspherical, "khvorostyanov_nonspherical"}
    };
  };
};
