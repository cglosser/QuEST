#ifndef LAGRANGE_SET_H
#define LAGRANGE_SET_H

#include <boost/multi_array.hpp>

namespace Interpolation {
  class UniformLagrangeSet;
  typedef boost::multi_array<double, 2> InterpolationTable;
}

class Interpolation::UniformLagrangeSet {
 public:

  const int order;
  InterpolationTable weights;

  UniformLagrangeSet(const int);
  UniformLagrangeSet(const double, const int);

  void calculate_weights(const double);
};

#endif
