// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/detail/bisect.hpp>
#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/henry.hpp>
#include <libcloudph++/common/dissoc.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct chem_volfun
      {
        const real_t pi;

        // ctor (pi() is not a __device__ function...)
        chem_volfun() :
          pi(boost::math::constants::pi<real_t>())
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2) const
        {
#if !defined(__NVCC__)
	  using std::pow;
#endif
          return real_t(4./3) * pi * (pow(rw2, real_t(3./2)));
        }
      };

      template <typename real_t>
      struct chem_absfun
      {
        const quantity<common::amount_over_volume_over_pressure, real_t> H;
        const quantity<common::mass_over_amount, real_t> M;
        const quantity<si::dimensionless, real_t> c;

        // ctor
        chem_absfun(
          const quantity<common::amount_over_volume_over_pressure, real_t> &H,
          const quantity<common::mass_over_amount, real_t> &M,
          const quantity<si::dimensionless, real_t> &c
        ) : 
          H(H), M(M), c(c)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &V, const real_t &p) const
        { 
          return (
            H                       // Henry's law
            * c * (p * si::pascals) // volume concentration -> partial pressure
            * V * si::cubic_metres  // drop volume 
	    * M                     // moles -> kilograms
          ) / si::kilograms;
        }
      };

      template <typename real_t>
      struct chem_minfun
      {
	const quantity<si::mass, real_t> &m_SO2, &m_S_VI;
	const quantity<si::volume, real_t> &V;
        
        // ctor
        BOOST_GPU_ENABLED
        chem_minfun(
	  const quantity<si::mass, real_t> &m_SO2,
	  const quantity<si::mass, real_t> &m_S_VI,
	  const quantity<si::volume, real_t> &V
        ) :
          m_SO2(m_SO2), m_S_VI(m_S_VI), V(V)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &arg) const
        {
	  using namespace common::molar_mass;
	  using namespace common::dissoc;

          const quantity<si::mass, real_t> m_H = arg * si::kilograms;

          return (-m_H + M_H<real_t>() * (
            // HSO3 to SO3 dissoctation 
            real_t(2) * // "2-" ion
            K_SO2<real_t>() * K_HSO3<real_t>() / M_SO2<real_t>() * m_SO2 
            * (M_H<real_t>() * M_H<real_t>()) 
            * (V * V) 
            / (m_H * m_H)
            +
            // H20*SO2 to HSO3 dissociation
            K_SO2<real_t>() / M_SO2<real_t>() * m_SO2 * M_H<real_t>() * V / m_H
            +
            // dissociation of pure water 
            K_H2O<real_t>() * M_H<real_t>() * (V*V) / m_H
            +
            // dissociation of S_VI to HSO4
            m_H * m_S_VI / V / M_H2SO4<real_t>() / M_H<real_t>()
	    / ((m_H) / M_H<real_t>() / V + K_HSO4<real_t>())
            +
            // dissociation of HSO4 to SO4
            real_t(2) * // "2-" ion
            K_HSO4<real_t>() * m_S_VI 
            / M_H2SO4<real_t>()
            / ((m_H) / M_H<real_t>() / V + K_HSO4<real_t>())


          )) / si::kilograms;
        }
      };

      template <typename real_t>
      struct chem_curie
      {
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::mass, real_t>
            m_H    = thrust::get<0>(tpl) * si::kilograms,
            m_SO2  = thrust::get<1>(tpl) * si::kilograms,
            m_S_VI = thrust::get<2>(tpl) * si::kilograms; 
          const quantity<si::volume, real_t> 
            V      = thrust::get<3>(tpl) * si::cubic_metres;

	  real_t tol = 44; // TODO!
          return common::detail::bisect(
	    detail::chem_minfun<real_t>(m_SO2, m_S_VI, V),
	    real_t(22), // min -> TODO! (pure water)
	    real_t(44), // max -> TODO
	    tol
	  ); 
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
      // 0/4: calculating drop volumes
      thrust_device::vector<real_t> &V(tmp_device_real_part);
      thrust::transform(
        rw2.begin(), rw2.end(),         // input
        V.begin(),                      // output 
        detail::chem_volfun<real_t>()   // op
      );

      // 1/4: equilibrium stuff: gas absortption
      // TODO: open/close system logic
      // TODO: K=K(T)
      using namespace common::henry;
      using namespace common::molar_mass;
      assert(SO2 == 0 && H2O2 == 1 && O3 == 2);
      static const quantity<common::amount_over_volume_over_pressure, real_t> H_[chem_gas_n] = {
	H_SO2<real_t>(), H_H2O2<real_t>(), H_O3<real_t>()
      };
      static const quantity<common::mass_over_amount, real_t> M_[chem_gas_n] = {
	M_SO2<real_t>(), M_H2O2<real_t>(), M_O3<real_t>()
      };
      for (int i = 0; i < chem_gas_n; ++i)
      {
	thrust::transform(
	  V.begin(), V.end(),                                        // input - 1st arg
	  thrust::make_permutation_iterator(p.begin(), ijk.begin()), // input - 2nd arg
	  che[i].begin(),                                            // output
	  detail::chem_absfun<real_t>(H_[i], M_[i], chem_gas[i])       // op
	);
      }

      // 2/4: equilibrium stuff: dissociation
      {
        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::iterator,
            typename thrust_device::vector<real_t>::iterator,
            typename thrust_device::vector<real_t>::iterator,
            typename thrust_device::vector<real_t>::iterator
          >
        > zip_it_t;

	thrust::transform(
	  zip_it_t(thrust::make_tuple(che[H].begin(), che[SO2].begin(), che[S_VI].begin(), V.begin())),  // input - begin
	  zip_it_t(thrust::make_tuple(che[H].end(),   che[SO2].end(),   che[S_VI].end(),   V.end())  ),  // input - end
	  che[H].begin(),                 // output
	  detail::chem_curie<real_t>()    // op
	);
      }

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
