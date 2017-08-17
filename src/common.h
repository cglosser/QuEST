#ifndef COMMON_H
#define COMMON_H

#include <boost/multi_array.hpp>
#include <complex>

typedef std::complex<double> cmplx;
constexpr cmplx iu(0, 1);

template <class T>
using Array = boost::multi_array<T, 2>;

template <class T>
using SpacetimeArray = boost::multi_array<T, 4>;

#endif
