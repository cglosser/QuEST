#ifndef COMMON_H
#define COMMON_H

#include <boost/multi_array.hpp>
#include <complex>

typedef std::complex<double> cmplx;
constexpr cmplx iu(0, 1);
typedef boost::multi_array<double, 2> DblArray;
typedef boost::multi_array<cmplx, 2> CmplxArray;

#endif
