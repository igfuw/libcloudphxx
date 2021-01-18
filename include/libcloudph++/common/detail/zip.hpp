/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Zip iterator implementation using C++11 variadic templates
  */

#pragma once

// http://stackoverflow.com/questions/8511035/sequence-zip-function-for-c11

#include <boost/iterator/zip_iterator.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/range.hpp>

template <typename... T>
auto zip(T&... containers)
  -> boost::iterator_range<
    boost::zip_iterator<
      decltype(
        std::make_tuple(std::begin(containers)...)
      )
    >
  >
{
  return boost::make_iterator_range(
    boost::make_zip_iterator(
      std::make_tuple(std::begin(containers)...)
    ),
    boost::make_zip_iterator(
      std::make_tuple(std::end(containers)...)
    )
  );
}
