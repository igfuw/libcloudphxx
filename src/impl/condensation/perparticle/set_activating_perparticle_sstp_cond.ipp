namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      struct if_activates
      {
        template<class tpl_t>
        BOOST_GPU_ENABLED bool operator()(tpl_t rw2_drw2_sstp_rc2) noexcept
        {
          const auto &rw2 = cuda::std::get<0>(rw2_drw2_sstp_rc2);
          const auto &drw2 = cuda::std::get<1>(rw2_drw2_sstp_rc2);
          const auto &sstp_cond = cuda::std::get<2>(rw2_drw2_sstp_rc2);
          const auto &rc2 = cuda::std::get<3>(rw2_drw2_sstp_rc2);
          // activation condition:
          if ( rw2 < rc2 && (rw2 + sstp_cond * drw2) > rc2 )
            return true;
          // deactivation condition:
          else if ( rw2 > rc2 && (rw2 + sstp_cond * drw2) < rc2 )
            return true;
          else
            return false;
        }
      };
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::set_activating_perparticle_sstp_cond(const unsigned int &n)
    {
      auto &perparticle_sstp_cond = perparticle_sstp_cond_gp->get();
      const auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();
      const auto &drw2 = drw2_gp->get();

      // computing rc2; TODO: same done in update_incloud_time and in particles_diag::diag_rw_ge_rc - refactor
      auto rc2_g = tmp_device_real_part.get_guard();
      thrust_device::vector<real_t> &rc2 = rc2_g.get();

      thrust::transform(
        rd3.begin(), rd3.end(), // input - 1st arg
        thrust::make_zip_iterator(make_tuple(
          kpa.begin(), 
          thrust::make_permutation_iterator(
            T.begin(),
            ijk.begin()
          )
        )),                                   // input - 2nd arg 
        rc2.begin(),                          // output
        detail::rw3_cr<real_t>()              // op
      );
      // thrust::fill(rc2.begin(), rc2.end(), real_t(1)); // test
      
      thrust::transform_if(
        thrust::make_constant_iterator<unsigned int>(n),
        thrust::make_constant_iterator<unsigned int>(n) + n_part,
        thrust::make_zip_iterator(
          thrust::make_tuple(
              rw2.begin(),
              drw2.begin(),
              perparticle_sstp_cond.begin(),
              rc2.begin()
          )
        ),
        perparticle_sstp_cond.begin(),
        cuda::std::identity(),
        detail::if_activates()
      );

      // thrust::copy_if(
      //   thrust::make_constant_iterator<unsigned int>(n),
      //   thrust::make_constant_iterator<unsigned int>(n) + n_part,
      //   unconverged_mask.begin(),
      //   perparticle_sstp_cond.begin(),
      //   cuda::std::identity()
      // );
    }
  };
};