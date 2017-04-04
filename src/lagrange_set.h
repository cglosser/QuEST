#ifndef LAGRANGE_SET_H
#define LAGRANGE_SET_H

#include <boost/multi_array.hpp>

namespace Interpolation {
  class UniformLagrangeSet;
  typedef boost::multi_array<double, 2> InterpolationTable;
  constexpr int NUM_DERIVATIVES = 3;
}

class Interpolation::UniformLagrangeSet {
 public:
  InterpolationTable weights;

  UniformLagrangeSet(const int);
  UniformLagrangeSet(const double, const int, const double = 1);

  void calculate_weights(const double, const double = 1);
  int order() const { return order_; }

 private:
  int order_;
};

#endif
