#pragma once
#include <thrust/tuple.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      //wrapper for the tuple tpl_rw_t needed by collision kernel
      template <class real_t, class n_t>
      struct tpl_calc_wrap
      {
        typedef thrust::tuple<
             n_t,           n_t,        // n   (multiplicity)
          real_t,        real_t,        // rw2 (wet radius squared)
          real_t,        real_t,        // vt  (terminal velocity)
          real_t,        real_t,        // rd3 (dry radius cubed)
          real_t,        real_t         // number of collisions
        > tpl_rw_t;

        typedef thrust::tuple<
          real_t,                       // rhod
          real_t                        // eta
        > tpl_ro_calc_t;

        tpl_rw_t tpl_rw;
        tpl_ro_calc_t tpl_ro_calc;

        BOOST_GPU_ENABLED
        tpl_calc_wrap(tpl_rw_t _tpl_rw, tpl_ro_calc_t _tpl_ro_calc):
          tpl_rw(_tpl_rw),
          tpl_ro_calc(_tpl_ro_calc)
          {}

        BOOST_GPU_ENABLED
        tpl_rw_t get_rw () const
        {
          return tpl_rw;
        }

        BOOST_GPU_ENABLED
        tpl_ro_calc_t get_ro_calc () const
        {
          return tpl_ro_calc;
        }
      };
    }
  }
}


