#ifndef LAGRANGE_SET_H
#define LAGRANGE_SET_H

#include <boost/multi_array.hpp>

class UniformLagrangeSet {
 public:
  typedef boost::multi_array<double, 2> Array;

  const int order;
  const double dt;
  Array weights;

  UniformLagrangeSet(const int, const double);
  UniformLagrangeSet(const double, const int, const double);

  void calculate_weights(const double);
};

#endif
