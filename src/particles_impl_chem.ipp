// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/detail/bisect.hpp>
#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/henry.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct chem_absfun
      {
        const quantity<common::amount_over_volume_over_pressure, real_t> H;
        const quantity<common::mass_over_amount, real_t> M;
        const quantity<si::dimensionless, real_t> c;
        const real_t pi;

        // ctor
        chem_absfun(
          const quantity<common::amount_over_volume_over_pressure, real_t> &H,
          const quantity<common::mass_over_amount, real_t> &M,
          const quantity<si::dimensionless, real_t> &c
        ) : 
          H(H), M(M), c(c), pi(boost::math::constants::pi<real_t>())  
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2, const real_t &p) const
        { 
#if !defined(__NVCC__)
	  using std::pow;
#endif
          return (
            H // Henry's law
            * c * (p * si::pascals) // volume concentration -> partial pressure
            * real_t(4./3) * pi * (pow(rw2, real_t(3./2))) * si::cubic_metres // drop volume
	    * M // moles -> kilograms
          ) / si::kilograms;
        }
      };

      template <typename real_t>
      struct chem_minfun
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &m_H) const
        {
          return 0;
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem(
      const real_t &dt,
      const std::vector<real_t> &chem_gas
    )
    {   
std::cerr << "@particles_t::impl::chem()" << std::endl;
      // 1/4: equilibrium stuff: gas absortption
      // TODO: open/close system logic
      // TODO: K=K(T)
      using namespace common::henry;
      using namespace common::molar_mass;
std::cerr << SO2 << std::endl;
std::cerr << H2O2 << std::endl;
std::cerr << O3 << std::endl;
      assert(SO2 == 0 && H2O2 == 1 && O3 == 2);
      static const quantity<common::amount_over_volume_over_pressure, real_t> H[chem_gas_n] = {
	H_SO2<real_t>(), H_H2O2<real_t>(), H_O3<real_t>()
      };
      static const quantity<common::mass_over_amount, real_t> M[chem_gas_n] = {
	M_SO2<real_t>(), M_H2O2<real_t>(), M_O3<real_t>()
      };
      for (int i = 0; i < chem_gas_n; ++i)
      {
	thrust::transform(
	  rw2.begin(), rw2.end(),                                    // input - 1st arg
	  thrust::make_permutation_iterator(p.begin(), ijk.begin()), // input - 2nd arg
	  che[i].begin(),                                            // output
	  detail::chem_absfun<real_t>(H[i], M[i], chem_gas[i])       // op
	);
      }

      // 2/4: equilibrium stuff: dissociation
      real_t tol = 44;
      real_t range = common::detail::bisect(
        detail::chem_minfun<real_t>(),
        real_t(22), // min -> pure water
        real_t(44), // max -> TODO
        tol
      ); 

      // 3/4: non-equilibrium stuff
      // TODO


      // 4/4: recomputing dry radii
      {
        using namespace thrust::placeholders;
        // TODO: using namespace for S_VI
        thrust::transform(
          che[S_VI].begin(), che[S_VI].end(),                      // input
          rd3.begin(),                                             // output
          (real_t(3./4) / pi<real_t>() / opts_init.chem_rho) * _1  // op
        );
      };
    }
  };  
};
