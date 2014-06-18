// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct chem_init
      {
        real_t mltpl;

        // ctor
        chem_init(const real_t &rho)
          : mltpl(real_t(4./3) * pi<real_t>() * rho)
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
std::cerr << "@init_chem()" << std::endl;
      // TODO: don't do it if not using chem...

      for (int i = chem_aq_begin; i < chem_aq_end; ++i)
      {
        // memory allocation
        che[i].resize(n_part);

std::vector<real_t> che_init(chem_aq_end+1);

	// initialising values 
	thrust::transform(
	  rw2.begin(), rw2.end(),        // input
	  che[i].begin(),                // output
	  detail::chem_init<real_t>(che_init[i]) // op
	);
      }
    }
  };
};
