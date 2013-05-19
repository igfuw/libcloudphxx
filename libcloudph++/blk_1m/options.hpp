#pragma once

namespace libcloudphxx
{
  namespace blk_1m
  {
    template<typename real_t>
    struct opts
    {
      bool 
        cevp = true, 
        revp = true, 
        conv = true, 
        clct = true, 
        sedi = true;
      real_t dt = 0;
    };
  }
};
