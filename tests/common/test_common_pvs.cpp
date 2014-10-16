//<listing-1>
#include <libcloudph++/common/const_cp.hpp>

int main()
{
  using namespace libcloudphxx::common::const_cp;
  assert(p_vs(273.16 * si::kelvins) == 611.73 * si::pascals);
}
//</listing-1>
