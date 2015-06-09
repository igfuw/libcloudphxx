#pragma once

#include <list>
// TODO: other includes

template <typename real_t>
using outmom_t = std::list<
  std::pair<
    std::pair<
      quantity<si::length, real_t>,
      quantity<si::length, real_t>
    >, 
    std::vector<int> 
  >
>;

// TODO: option parsing code should go here...
