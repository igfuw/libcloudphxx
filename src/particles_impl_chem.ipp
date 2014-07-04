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
            // H2O*SO2 to HSO3 dissociation
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
      struct chem_curie_pH // TODO: does it have to be a struct/functor - perhaps ordinary function would suffice?
      {
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) const
        {
          using namespace common::molar_mass;

          const quantity<si::mass, real_t>
            m_SO2  = thrust::get<0>(tpl) * si::kilograms,
            m_S_VI = thrust::get<1>(tpl) * si::kilograms; 
          const quantity<si::volume, real_t> 
            V      = thrust::get<2>(tpl) * si::cubic_metres;

	  real_t tol = 1e-44; // TODO!
          real_t m_H_pure = ((real_t(1e-7 * 1e3) * si::moles / si::cubic_metres) * V * M_H<real_t>()) / si::kilograms;
          real_t m_H = common::detail::bisect(
	    detail::chem_minfun<real_t>(m_SO2, m_S_VI, V),
            m_H_pure, // min -> (pure water)
	    real_t(1e-10), // max -> TODO
	    tol
	  ); 
          //std::cerr << "  " << m_H_pure << " ... " << m_H << std::endl;
          // TODO: asserts for K = f(m_H, m_...)
          return m_H;
        }
      };
 
      template <typename real_t, int plnk>
      struct chem_curie_diag // TODO: does it have to be a struct/functor - perhaps ordinary function would suffice?
      {
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::volume, real_t> 
            V      = thrust::get<0>(tpl) * si::cubic_metres;
          const quantity<si::mass, real_t>
            m_H    = thrust::get<1>(tpl) * si::kilograms,
            m_SO2  = thrust::get<2>(tpl) * si::kilograms,
            m_HSO3 = thrust::get<3>(tpl) * si::kilograms,
            m_S_VI = thrust::get<4>(tpl) * si::kilograms;

          switch (plnk)
          {
	    using namespace common::dissoc;     // K-prefixed
	    using namespace common::molar_mass; // M-prefixed

            case OH:
              return (
                M_H<real_t>() * M_OH<real_t>() 
                * K_H2O<real_t>() // note: dissociation constant for pure water is actually k*[H2O] (Seinfeld and Pandis p 345)
                * V * V / m_H
              ) / si::kilograms; 
            case HSO3:
              return (
                V / m_H * m_SO2
                * M_HSO3<real_t>() * K_SO2<real_t>() * M_H<real_t>() / M_SO2<real_t>()
              ) / si::kilograms;
            case SO3:
              return (
                V
                * M_SO3<real_t>()
                * K_HSO3<real_t>() 
                * m_HSO3 / m_H
                * M_H<real_t>() / M_HSO3<real_t>()
              ) / si::kilograms;
            case HSO4:
              return (
                M_HSO4<real_t>() / M_H<real_t>() / M_H2SO4<real_t>()
                * m_H
                * m_S_VI
                / V 
                / ( m_H / M_H<real_t>() / V + K_HSO4<real_t>())
              ) / si::kilograms;
            case SO4:
              return (
                M_SO4<real_t>() / M_H2SO4<real_t>()
                * K_HSO4<real_t>() 
                * m_S_VI
                / (m_H / M_H<real_t>() / V + K_HSO4<real_t>())
              ) / si::kilograms;
            default:
              assert(false);
          }
        }
      };

      template <typename real_t>
      BOOST_GPU_ENABLED
      void chem_rhs(
        const thrust_device::vector<real_t> &psi, 
        thrust_device::vector<real_t> &dot_psi,
        const real_t &dt_
      )
      {
        const quantity<si::time, real_t> dt = dt_;
        
      }
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem(
      const real_t &dt,
      const std::vector<real_t> &chem_gas
    )
    {   
      using namespace common::henry;      // H-prefixed
      using namespace common::molar_mass; // M-prefixed

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
	  chem_bgn[i],                                               // output
	  detail::chem_absfun<real_t>(H_[i], M_[i], chem_gas[i])     // op
	);
      }

      // 2/4: equilibrium stuff: dissociation
      // H+ 
      {
        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::iterator, // SO2
            typename thrust_device::vector<real_t>::iterator, // S_VI
            typename thrust_device::vector<real_t>::iterator  // V
          >
        > zip_it_t;

	thrust::transform(
	  zip_it_t(thrust::make_tuple(chem_bgn[SO2], chem_bgn[S_VI], V.begin())),  // input - begin
	  zip_it_t(thrust::make_tuple(chem_end[SO2], chem_end[S_VI], V.end())  ),  // input - end
	  chem_bgn[H],                                                             // output
	  detail::chem_curie_pH<real_t>()                                          // op
	);
      }
      {
        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::iterator, // V
            typename thrust_device::vector<real_t>::iterator, // H 
            typename thrust_device::vector<real_t>::iterator, // SO2
            typename thrust_device::vector<real_t>::iterator, // HSO3
            typename thrust_device::vector<real_t>::iterator  // S_VI
          >
        > zip_it_t;

        zip_it_t 
          arg_begin(thrust::make_tuple(V.begin(), chem_bgn[H], chem_bgn[SO2], chem_bgn[HSO3], chem_bgn[S_VI])),
          arg_end(  thrust::make_tuple(V.end(),   chem_end[H], chem_end[SO2], chem_end[HSO3], chem_end[S_VI]));
        
        thrust::transform(arg_begin, arg_end, chem_bgn[OH  ], detail::chem_curie_diag<real_t, OH  >());
        thrust::transform(arg_begin, arg_end, chem_bgn[HSO3], detail::chem_curie_diag<real_t, HSO3>()); // note: has to be computed before SO3
        thrust::transform(arg_begin, arg_end, chem_bgn[SO3 ], detail::chem_curie_diag<real_t, SO3 >());
        thrust::transform(arg_begin, arg_end, chem_bgn[HSO4], detail::chem_curie_diag<real_t, HSO4>());
        thrust::transform(arg_begin, arg_end, chem_bgn[SO4 ], detail::chem_curie_diag<real_t, SO4 >());
      }

      // 3/4: non-equilibrium stuff
      {
/*
        chem_stepper.do_step(
          detail::chem_rhs<real_t>, 
          chem_state, 
          0,
          dt
        );
*/
      }


      // 4/4: recomputing dry radii
      {
        namespace arg = thrust::placeholders;
        // TODO: using namespace for S_VI
        thrust::transform(
          chem_bgn[S_VI], chem_end[S_VI],                               // input
          rd3.begin(),                                                  // output
          (real_t(3./4) / pi<real_t>() / opts_init.chem_rho) * arg::_1  // op
        );
      };
    }
  };  
};
