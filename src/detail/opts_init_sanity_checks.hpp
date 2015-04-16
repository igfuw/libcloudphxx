
#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class real_t>
      void opts_init_sanity_checks(const opts_init_t<real_t> &opts_init)
      {
        if(opts_init.terminal_velocity == undefined) throw std::runtime_error("Please define terminal velocity formula to use");
      }
    };
  };
};
