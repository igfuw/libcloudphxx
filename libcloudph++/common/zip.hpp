#pragma once

// http://stackoverflow.com/questions/8511035/sequence-zip-function-for-c11

#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

template <typename... T>
auto zip(T&... containers) 
  -> boost::iterator_range<
    boost::zip_iterator<
      decltype(
        boost::make_tuple(std::begin(containers)...)
      )
    >
  >
{
  return boost::make_iterator_range(
    boost::make_zip_iterator(
      boost::make_tuple(std::begin(containers)...)
    ),
    boost::make_zip_iterator(
      boost::make_tuple(std::end(containers)...)
    )
  );
}
