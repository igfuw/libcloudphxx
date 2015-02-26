#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      //wrapper for the tuple tpl_rw_t needed by collision kernel
      template <class real_t, class n_t>
      struct tpl_rw_t_wrap
      {
        typedef thrust::tuple<
             n_t,           n_t,           // n   (multiplicity)
          real_t,        real_t,        // rw2 (wet radius squared)
          real_t,        real_t,        // vt  (terminal velocity)
          real_t,        real_t         // rd3 (dry radius cubed)
        > tpl_rw_t;

        tpl_rw_t tpl_rw;

        BOOST_GPU_ENABLED
        tpl_rw_t_wrap(const tpl_rw_t &_tpl_rw):
          tpl_rw(_tpl_rw)
          {}

        BOOST_GPU_ENABLED
        tpl_rw_t operator ()() const
        {
          return tpl_rw;
        }
      };
    }
  }
}


