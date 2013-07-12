// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      namespace detail
      {

        template <typename real_t>
	struct eval_and_multiply
	{ 
          const unary_function<real_t> &fun;
          const real_t &mul;

          // ctor
          eval_and_multiply(
            const unary_function<real_t> &fun, 
            const real_t &mul
          ) 
            : fun(fun), mul(mul)
          {}

	  real_t operator()(real_t x) 
	  { 
	    return mul * fun(x); 
	  } 
        };

      };
    };
  };
};
