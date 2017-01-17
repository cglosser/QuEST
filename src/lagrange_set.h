#ifndef LAGRANGE_SET_H
#define LAGRANGE_SET_H

#include <boost/multi_array.hpp>

namespace Interpolation {
  class UniformLagrangeSet;
}

class Interpolation::UniformLagrangeSet {
 public:
  typedef boost::multi_array<double, 2> Array;

  const int order;
  Array weights;

  UniformLagrangeSet(const int);
  UniformLagrangeSet(const double, const int);

  void calculate_weights(const double);
};

#endif
