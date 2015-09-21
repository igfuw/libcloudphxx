#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      template<class real_t>
      int get_dev_nx(const opts_init_t<real_t> &opts_init, const int &dev_no)
      {
        if(dev_no < opts_init.dev_count-1)
          return opts_init.nx / opts_init.dev_count + .5;
        else
          return opts_init.nx - dev_no * int(opts_init.nx / opts_init.dev_count + .5); 
      } 
    }
  }
}
