// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/molar_mass.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct chem_init_water // water
      {
        const real_t mltpl;

        // ctor
        chem_init_water(const real_t &pH, const quantity<common::mass_over_amount, real_t> &M) : 
          mltpl(pow(10,-pH) * real_t(4./3) * (M / si::kilograms * si::moles) * pi<real_t>())
        {}
 
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2)
        {
          return mltpl * pow(rw2, real_t(3./2));
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_chem()
    {
      // TODO: don't do it if not using chem...

      // memory allocation
      chem_bgn.resize(chem_aq_n);
      chem_end.resize(chem_aq_n);
      chem_noneq.resize(chem_rhs_n * n_part);
      chem_equil.resize((chem_aq_n - chem_rhs_n) * n_part);
      chem_stepper.adjust_size(chem_noneq);
      
      // helper iterators
      for (int i = 0; i < chem_aq_n; ++i)
      {
        thrust_device::vector<real_t> &vec(
          i < chem_rhs_n
            ? chem_noneq
            : chem_equil
        );
        const int off = 
          i < chem_rhs_n
            ? 0
            : chem_rhs_n;
        chem_bgn[i] = vec.begin() + (i   - off) * n_part;
        chem_end[i] = vec.begin() + (i+1 - off) * n_part;
      }

      // initial values
      for (int i = 0; i < chem_aq_n; ++i)
      {
        switch (i)
        {
          case OH:
          case H:
            // pH = 7 
            {
              quantity<common::mass_over_amount, real_t> M;
              switch (i)
              {
                case OH: M = common::molar_mass::M_OH<real_t>(); break;
                case  H: M = common::molar_mass::M_H<real_t>(); break;
                default: assert(false);
              }
	      thrust::transform(
		rw2.begin(), rw2.end(),               // input
		chem_bgn[i],                          // output
		detail::chem_init_water<real_t>(7, M) // op
	      );
            }
            break;
          case S_VI:
            // implied by the lognormal distro
            {
              namespace arg = thrust::placeholders;
	      thrust::transform(
		rd3.begin(), rd3.end(),                                      // input
		chem_bgn[i],                                                 // output
		(real_t(4./3) * pi<real_t>() * opts_init.chem_rho) * arg::_1 // op
	      );
            }
            break;
          default: 
            // ... TODO: epsilon?
            thrust::fill(chem_bgn[i], chem_end[i], 0);
        }
      }
    }
  };
};
