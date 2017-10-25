// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */


// contains definitions of members of particles_t specialized for multiple GPUs
#include <future>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <class real_t>
      std::map<output_t, real_t> empty_out_map()
      {
        std::map<output_t, real_t> res;
        for(int i=0; i < chem_all+2; ++i) 
          res[static_cast<output_t>(i)] = 0.;
        return res;
      }
    }
    // diagnostic methods
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_RH()
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_RH);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_vel_div()
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_vel_div);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_sd_conc()
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_sd_conc);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_dry_rng, r_mi, r_mx);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_wet_rng, r_mi, r_mx);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_mom(const int &k)
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_dry_mom, k);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mom(const int &k)
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_wet_mom, k);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mass_dens(const real_t &a, const real_t &b)
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_wet_mass_dens, a, b);
    }

    // ...
//</listing>

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_chem(const enum chem_species_t &spec)
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_chem, spec);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_rw_ge_rc()
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_rw_ge_rc);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_RH_ge_Sc()
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_RH_ge_Sc);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_all()
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_all);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_precip_rate()
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_precip_rate);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_max_rw()
    {
      mcuda_run(&particles_t<real_t, CUDA>::diag_max_rw);
    }

    template <typename real_t>
    real_t* particles_t<real_t, multi_CUDA>::outbuf()
    {
      // run fill_outbuf on each gpu
      std::vector<std::thread> threads;
      for (int i = 0; i < glob_opts_init.dev_count; ++i)
      {
        threads.emplace_back(
          detail::set_device_and_run, i, 
          std::bind(
            &particles_t<real_t, CUDA>::impl::fill_outbuf,
            &(*(particles[i]->pimpl))
          )
        );
      }
      for (auto &th : threads) th.join();

      for(auto &p : particles) // TODO: perform this copy in parallell?
      {
        thrust::copy(
          p->pimpl->tmp_host_real_cell.begin(),
          p->pimpl->tmp_host_real_cell.end(),
          real_n_cell_tot.begin() + p->pimpl->n_cell_bfr
        );
      }
      return &(*(real_n_cell_tot.begin()));
    }

    template<class real_t>
    std::map<output_t, real_t> add_puddle(std::map<output_t, real_t> x, std::map<output_t, real_t> y){
      std::map<output_t, real_t> res;
      for(int i=0; i < chem_all+2; ++i) 
        res[static_cast<output_t>(i)] = x[static_cast<output_t>(i)] + y[static_cast<output_t>(i)];
      return res;
    }

    template <typename real_t>
    std::map<output_t, real_t> particles_t<real_t, multi_CUDA>::diag_puddle()
    {
      using pudmap_t = std::map<output_t, real_t>;
      pudmap_t res = detail::empty_out_map<real_t>();

      std::vector<std::future<pudmap_t>> futures(glob_opts_init.dev_count);
      for (int i = 0; i < glob_opts_init.dev_count; ++i)
      {
        futures[i] = std::async(
          std::launch::async,
          [i, this](){
            gpuErrchk(cudaSetDevice(i));
            return this->particles[i]->diag_puddle();
          }
        );
      }
      // TODO: optimize this...
      for (int i = 0; i < glob_opts_init.dev_count; ++i)
      {
        res = add_puddle(res, futures[i].get());
      }
      return res;
    }
  };
};
