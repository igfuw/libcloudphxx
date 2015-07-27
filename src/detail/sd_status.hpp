#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      enum sd_stat_t { active, inactive, to_init, to_rcyc }; 

      struct active_or_to_init : std::unary_function<bool, sd_stat_t>
      {
        BOOST_GPU_ENABLED
        bool operator()(const sd_stat_t &stat)
        {
          return (stat == active || stat == to_init);
        }
      };
    };
  };
};
