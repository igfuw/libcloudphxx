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
      std::map<common::output_t, real_t> empty_out_map()
      {
        std::map<common::output_t, real_t> res;
        for(int i=0; i < chem_all+2; ++i) 
          res[static_cast<common::output_t>(i)] = 0.;
        return res;
      }
    }
    // diagnostic methods
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_pressure()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_pressure);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_temperature()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_temperature);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_RH()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_RH);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_vel_div()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_vel_div);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_sd_conc()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_sd_conc);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_dry_rng, r_mi, r_mx);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_wet_rng, r_mi, r_mx);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_kappa_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_kappa_rng, r_mi, r_mx);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_ice()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_ice);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_water()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_water);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_rng_cons(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_dry_rng_cons, r_mi, r_mx);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_rng_cons(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_wet_rng_cons, r_mi, r_mx);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_kappa_rng_cons(
      const real_t &r_mi, const real_t &r_mx
    )
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_kappa_rng_cons, r_mi, r_mx);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_ice_cons()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_ice_cons);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_water_cons()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_water_cons);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_mom(const int &k)
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_dry_mom, k);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mom(const int &k)
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_wet_mom, k);
    }
 
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_kappa_mom(const int &k)
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_kappa_mom, k);
    }   
 
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_up_mom(const int &k)
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_up_mom, k);
    }   
 
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_vp_mom(const int &k)
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_vp_mom, k);
    }   
 
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wp_mom(const int &k)
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_wp_mom, k);
    }   
 
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_incloud_time_mom(const int &k)
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_incloud_time_mom, k);
    }   

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mass_dens(const real_t &a, const real_t &b)
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_wet_mass_dens, a, b);
    }

    // ...
//</listing>

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_chem(const enum chem_species_t &spec)
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_chem, spec);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_rw_ge_rc()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_rw_ge_rc);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_RH_ge_Sc()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_RH_ge_Sc);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_all()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_all);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_precip_rate()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_precip_rate);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_max_rw()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_max_rw);
    }

    template <typename real_t>
    real_t* particles_t<real_t, multi_CUDA>::outbuf()
    {
      // run fill_outbuf on each gpu
      std::vector<std::thread> threads;
      for (int i = 0; i < this->opts_init->dev_count; ++i)
      {
        auto &ipimpl = pimpl->particles[i]->pimpl;
        reset_guardp(ipimpl->outbuf_gp, ipimpl->tmp_host_real_cell);
        thrust::host_vector<real_t> &outbuf(ipimpl->outbuf_gp->get());

        threads.emplace_back(
          detail::set_device_and_run, i, 
          std::bind(
            &particles_t<real_t, CUDA>::impl::fill_outbuf,
            &(*(pimpl->particles[i]->pimpl)),
            outbuf
          )
        );
      }
      for (auto &th : threads) th.join();

      for(auto &p : pimpl->particles) // TODO: perform this copy in parallell?
      {
        thrust::host_vector<real_t> &outbuf(p->pimpl->outbuf_gp->get());

        thrust::copy(
          outbuf.begin(),
          outbuf.end(),
          pimpl->real_n_cell_tot.begin() + p->pimpl->n_cell_bfr
        );
        p->pimpl->outbuf_gp.reset(); // release the guard
      }
      return &(*(pimpl->real_n_cell_tot.begin()));
    }

    template<class real_t>
    std::map<common::output_t, real_t> add_puddle(std::map<common::output_t, real_t> x, std::map<common::output_t, real_t> y){
      std::map<common::output_t, real_t> res;
      for(int i=0; i < common::output_names.size(); ++i) 
        res[static_cast<common::output_t>(i)] = x[static_cast<common::output_t>(i)] + y[static_cast<common::output_t>(i)];
      return res;
    }

    template <typename real_t>
    std::map<common::output_t, real_t> particles_t<real_t, multi_CUDA>::diag_puddle()
    {
      using pudmap_t = std::map<common::output_t, real_t>;
      pudmap_t res = detail::empty_out_map<real_t>();

      std::vector<std::future<pudmap_t>> futures(this->opts_init->dev_count);
      for (int i = 0; i < this->opts_init->dev_count; ++i)
      {
        futures[i] = std::async(
          std::launch::async,
          [i, this](){
            gpuErrchk(cudaSetDevice(i));
            return this->pimpl->particles[i]->diag_puddle();
          }
        );
      }
      // TODO: optimize this...
      for (int i = 0; i < this->opts_init->dev_count; ++i)
      {
        res = add_puddle(res, futures[i].get());
      }
      return res;
    }


    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_revp20()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_revp20);
    }
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_revp25()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_revp25);
    }
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_revp32()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_revp32);
    }
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_accr20()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_accr20);
    }
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_accr25()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_accr25);
    }
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_accr32()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_accr32);
    }
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_acnv20()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_acnv20);
    }
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_acnv25()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_acnv25);
    }
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_acnv32()
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::diag_acnv32);
    }
  };
};
