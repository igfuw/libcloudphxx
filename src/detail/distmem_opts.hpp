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

      // adjust opts_int for a distributed memory system
      // returns n_x_bfr
      template <class real_t>
      int distmem_opts(opts_init_t<real_t> &opts_init, const int &rank, const int &size)
      {
        int n_x_bfr = rank * get_dev_nx(opts_init, 0);

        opts_init.nx = detail::get_dev_nx(opts_init, rank);
 
        if(rank != 0)      opts_init.x0 = 0.;  // TODO: what if x0 greater than domain of first device?
        if(rank != size-1) opts_init.x1 = opts_init.nx * opts_init.dx;
        else               opts_init.x1 = opts_init.x1 - n_x_bfr * opts_init.dx;

        opts_init.n_sd_max = opts_init.n_sd_max / size + 1;

        return n_x_bfr;
      }
    }
  }
}
