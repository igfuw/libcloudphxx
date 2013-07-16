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
      struct exp3x
      { 
	__device__ 
	real_t operator()(real_t x) 
	{ 
	  return exp(3*x); 
	} 
      };


      template <typename T>
      struct equals_zero
      { 
	__device__ 
	bool operator()(T x) 
	{ 
	  return x == 0; 
	} 
      };

    };
  };
};
