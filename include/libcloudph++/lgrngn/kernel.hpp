#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
//<listing>
    enum class kernel_t { undefined, geometric, golovin, hall, hall_davis_no_waals, Long, onishi_hall, onishi_hall_davis_no_waals, hall_pinsky_1000mb_grav, hall_pinsky_cumulonimbus, hall_pinsky_stratocumulus, vohl_davis_no_waals}; 
//</listing>

    const std::unordered_map<kernel_t, std::string> kernel_name = {
      {kernel_t::undefined, "undefined"},
      {kernel_t::geometric, "geometric"},
      {kernel_t::golovin, "golovin"},
      {kernel_t::hall, "hall"},
      {kernel_t::hall_davis_no_waals, "hall_davis_no_waals"},
      {kernel_t::Long, "Long"},
      {kernel_t::onishi_hall, "onishi_hall"},
      {kernel_t::onishi_hall_davis_no_waals, "onishi_hall_davis_no_waals"},
      {kernel_t::hall_pinsky_1000mb_grav, "hall_pinsky_1000mb_grav"},
      {kernel_t::hall_pinsky_cumulonimbus, "hall_pinsky_cumulonimbus"},
      {kernel_t::hall_pinsky_stratocumulus, "hall_pinsky_stratocumulus"},
      {kernel_t::vohl_davis_no_waals, "vohl_davis_no_waals"}
    };
  };
};
