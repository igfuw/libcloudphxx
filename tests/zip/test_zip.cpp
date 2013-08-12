#include <vector>
#include <tuple>
#include <libcloudph++/common/detail/zip.hpp>

int main()
{
  std::vector<int> v;

  // for 1-10 elements works out-of-the box
  for (const auto &item : zip(v, v)) 
    ;

  // for >10 elementes it's trickier :(
  for (const auto &item : zip(v, v, v, v, v, v, v))
    ;
}
